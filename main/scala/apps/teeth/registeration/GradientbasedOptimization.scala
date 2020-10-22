package apps.teeth.registeration

import java.io.File

import apps.teeth.utilities.Paths.generalPath
import apps.util.AlignmentTransforms
import breeze.linalg.DenseVector
import scalismo.common.NearestNeighborInterpolator
import scalismo.geometry.{EuclideanVector, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, LBFGSOptimizer}
import scalismo.registration.{GaussianProcessTransformationSpace, L2Regularizer, MeanSquaresMetric, Registration}
import scalismo.statisticalmodel.LowRankGaussianProcess
import scalismo.ui.api.ScalismoUI

object GradientbasedOptimization {
  implicit val rng = scalismo.utils.Random(42)
  scalismo.initialize()

  val ui = ScalismoUI()

  val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath +"/Ref/model/tooth_gp_model_50000-components.h5")).get
  val modelGroup = ui.createGroup("model")
  val modelLmsFile = new File(generalPath, "/Ref/landmarks/ref_landmarks_7.json")
  val modelLms = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get
  val gpView = ui.addTransformation(modelGroup, model.gp, "gp")

  def main(args: Array[String]):Unit ={
     val sampleMeshesFileList = new java.io.File(generalPath+"/aligned/meshes/toBeRegistered").listFiles.sortBy(f => f.getName())
  //  val registeredMeshesFileList = new java.io.File(generalPath + "/Registered/mesh/bad_reg/").listFiles

    val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
    val lowRankGP = model.gp.interpolate(interpolator)
    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)

    ui.show(modelGroup,model, "moving mesh")


    val initialCoefficients = DenseVector.zeros[Double](model.rank)
    val registrationParameters = Seq(

      RegistrationParameters(regularizationWeight = 0.5, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.1, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.05, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.01, numberOfIterations = 1000, numberOfSampledPoints = 1000),

      RegistrationParameters(regularizationWeight = 0.005, numberOfIterations = 1000, numberOfSampledPoints = 2000),
      RegistrationParameters(regularizationWeight = 0.001, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.0005, numberOfIterations = 1000, numberOfSampledPoints = 1000),

        RegistrationParameters(regularizationWeight = 0.0001, numberOfIterations = 1000, numberOfSampledPoints = 1000),

          RegistrationParameters(regularizationWeight = 0.00001, numberOfIterations = 1000, numberOfSampledPoints = 1000)

    )
    sampleMeshesFileList.map{ m =>
      val f = MeshIO.readMesh(m).get

    //  ui.show(f, s"Fixed mesh ${m.getName}")
      val targetLmsFile = new File(generalPath, s"aligned/landmarks/${m.getName.replace(".stl",".json")}")
      val targetLmsInit = LandmarkIO.readLandmarksJson[_3D](targetLmsFile).get

      val center: EuclideanVector[_3D] = model.mean.pointSet.points.map(_.toVector).reduce(_ + _) * 1.0 / model.mean.pointSet.numberOfPoints.toDouble

      val targetTransform = AlignmentTransforms.computeTransform( targetLmsInit,modelLms, center.toPoint)
      val fixedMesh = f.transform(targetTransform)
      ui.show(fixedMesh, s"Fixed mesh ${m.getName}")
      println(s"Fixed mesh ${m.getName}")
      val finalCoefficients: DenseVector[Double] = registrationParameters.foldLeft(initialCoefficients)((modelCoefficients, regParameters) =>
        doRegistration(lowRankGP, model.referenceMesh,fixedMesh, regParameters, modelCoefficients))

  /*    val registrationTransformation = transformationSpace.transformForParameters(finalCoefficients)
      val targetMeshOperations = fixedMesh
      val projection = (pt : Point[_3D]) => {
        targetMeshOperations.pointSet.findClosestPoint(pt).point
      }

      val finalTransformation = registrationTransformation.andThen(projection)
     */
     // val projectedMesh = model.referenceMesh.transform(finalTransformation)
      val projectedMesh = model.instance(finalCoefficients)
      val resultGroup = ui.createGroup("result")
      val projectionView = ui.show(resultGroup, projectedMesh, s"projection${m.getName}")
      MeshIO.writeVTK(projectedMesh, new File(generalPath +s"/Registered/mesh/${m.getName.replace(".stl",".vtk")}"))

    }

  }

  case class RegistrationParameters(regularizationWeight : Double, numberOfIterations : Int, numberOfSampledPoints : Int)
  def doRegistration(
                      lowRankGP : LowRankGaussianProcess[_3D, EuclideanVector[_3D]],
                      referenceMesh : TriangleMesh[_3D],
                      targetmesh : TriangleMesh[_3D],
                      registrationParameters : RegistrationParameters,
                      initialCoefficients : DenseVector[Double]
                    ) : DenseVector[Double] =
  {
    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)
    val fixedImage = referenceMesh.operations.toDistanceImage
    val movingImage = targetmesh.operations.toDistanceImage
    val sampler1 = FixedPointsUniformMeshSampler3D(
      referenceMesh,
      registrationParameters.numberOfSampledPoints
    )


    val metric1 = MeanSquaresMetric(
      fixedImage,
      movingImage,
      transformationSpace,
      sampler1
    )
    val optimizer1 = LBFGSOptimizer(registrationParameters.numberOfIterations)
    val regularizer1 = L2Regularizer(transformationSpace)
    val registration = Registration(
      metric1,
      regularizer1,
      registrationParameters.regularizationWeight,
      optimizer1
    )
    val registrationIterator = registration.iterator(initialCoefficients)
    val visualizingRegistrationIterator = for ((it, itnum) <- registrationIterator.zipWithIndex) yield {
      println(s"object value in iteration $itnum is ${it.value}")
      gpView.coefficients = it.parameters
      it
    }
    val registrationResult = visualizingRegistrationIterator.toSeq.last
    registrationResult.parameters
  }

}
