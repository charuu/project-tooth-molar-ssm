/*
 *  Copyright University of Basel, Graphics and Vision Research Group
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package apps.teeth.registeration

import java.awt.Color
import java.io.File

import api.other._
import api.sampling._
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.{LoadTestData, alignShapes}
import scalismo.geometry._3D
import scalismo.io.MeshIO
import scalismo.mesh.{TriangleMesh, TriangleMesh3D}
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.{ScalismoUI, StatisticalMeshModelViewControls}

object IcpRegistrationSurfacevsNormal {
  val meshType = "Partial"
  def fitting(model: StatisticalMeshModel, targetMesh: TriangleMesh3D, numOfSamplePoints: Int, numOfIterations: Int,regType:String, showModel: Option[StatisticalMeshModelViewControls], initialParameters: Option[ModelFittingParameters] = None): TriangleMesh[_3D] = {

    val initPars =
      if (initialParameters.isDefined) {
        Option(initialParameters.get.shapeParameters.parameters, ModelFittingParameters.poseTransform(initialParameters.get))
      }
      else {
        None
      }

    val best = if(regType == "Normal") {
      val icpFitting = IcpBasedNormalFitting(model, targetMesh, numOfSamplePoints, projectionDirection = ModelAndTargetSampling, showModel = showModel)
      icpFitting.runfitting(numOfIterations, iterationSeq = Seq(2.0,1.0,0.1,0.01,0.001), initialModelParameters = initPars)

    }else{
      val icpFitting = IcpBasedSurfaceFitting(model, targetMesh, numOfSamplePoints, projectionDirection = TargetSampling, showModel = showModel)
      icpFitting.runfitting(numOfIterations, iterationSeq = Seq(1.0,0.1,0.01,0.001), initialModelParameters = initPars)

    }

    val t0 = System.currentTimeMillis()

    val t1 = System.currentTimeMillis()
    println(s"ICP-Timing: ${(t1 - t0) / 1000.0} sec")
    best
  }


  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    println(s"Starting Standard NonRigid ICP registrations!")

    val (modelD, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget(meshType)
    val alignedMeshLms = alignShapes(modelLms,targetMeshes,targetLms).align()

    alignedMeshLms.map { targetMeshLms =>
      val model = modelD.decimate(1000)

    val targetMesh = targetMeshLms._1._1.operations.decimate(1000)
    val targetName = targetMeshLms._1._2
    val targetLms = targetMeshLms._2._1

    val numOfEvaluatorPoints = model.referenceMesh.pointSet.numberOfPoints // Used for the evaluation
    val numOfIterations = 75 // Number of ICP iterations

    val ui = ScalismoUI(s"ICP-registration")
    val modelGroup = ui.createGroup("modelGroup")
    val targetGroup = ui.createGroup("targetGroup")
    val finalGroup = ui.createGroup("finalGroup")


     /* val commonLmNames = modelLms.map(_.id) intersect targetMeshLms._2._1.map(_.id)
      val newcorr = commonLmNames.map(name => (model.mean.pointSet.findClosestPoint(modelLms.find(_.id == name).get.point).id, targetMeshLms._2._1.find(_.id == name).get.point, true))

      val regressionData = newcorr.map { correspondence =>
        val littleNoise = SurfaceNoiseHelpers.surfaceNormalDependantNoise(model.referenceMesh.vertexNormals.atPoint(correspondence._1), 0.1, 0.1)

        (correspondence._1, correspondence._2, littleNoise)
      }

      val posterior = model.posterior(regressionData.toIndexedSeq)

*/
      val showModel = ui.show(modelGroup, model, "model")
      ui.show(modelGroup, modelLms, "landmarks")
      val showTarget = ui.show(targetGroup, targetMesh, "target")
      ui.show(targetGroup, targetLms, "landmarks")
      showTarget.color = Color.YELLOW

 //   val bestSurfaceRegistration = fitting(posterior, targetMesh, numOfEvaluatorPoints, numOfIterations = numOfIterations, "Surface", showModel = Option(showModel))
 //     ui.show(finalGroup, bestSurfaceRegistration, s"best-surface-fit-${targetName}")
    val bestNormalRegistration= fitting(model, targetMesh, numOfEvaluatorPoints, numOfIterations = numOfIterations, "Normal", showModel = Option(showModel))
      ui.show(finalGroup, bestNormalRegistration, s"best-normal-fit-${targetName}")

    //  MeshIO.writeVTK(bestSurfaceRegistration, new File(generalPath + s"/Registered/partialMeshes/${targetName}_surface.vtk"))
      MeshIO.writeVTK(bestNormalRegistration, new File(generalPath + s"/Registered/partialMeshes/${targetName}_normal.vtk"))

  //  RegistrationComparison.evaluateReconstruction2GroundTruth("SAMPLE", bestSurfaceRegistration, targetMesh)
    RegistrationComparison.evaluateReconstruction2GroundTruth("SAMPLE", bestNormalRegistration, targetMesh)
  }
  }
}
