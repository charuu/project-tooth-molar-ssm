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

package api.sampling.proposals

import java.io.File

import api.other.{IcpProjectionDirection, ModelSampling, TargetSampling}
import api.sampling.{ModelFittingParameters, ShapeParameters, SurfaceNoiseHelpers}
import apps.teeth.utilities.Paths.generalPath
import scalismo.common.{Field, NearestNeighborInterpolator, PointId}
import scalismo.geometry._
import scalismo.io.LandmarkIO
import scalismo.mesh.{TriangleMesh, TriangleMesh3D}
import scalismo.numerics.{Sampler, UniformMeshSampler3D}
import scalismo.registration.RigidTransformation
import scalismo.sampling.{ProposalGenerator, TransitionProbability}
import scalismo.statisticalmodel.{LowRankGaussianProcess, MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.utils.{Memoize, Random}

import scala.collection.immutable.Stream.Empty

case class NonRigidIcpPartialProposal(
                                model: StatisticalMeshModel,
                                target: TriangleMesh3D,
                                stepLength: Double,
                                tangentialNoise: Double,
                                noiseAlongNormal: Double,
                                numOfSamplePoints: Int,
                                projectionDirection: IcpProjectionDirection = ModelSampling,
                                boundaryAware: Boolean = true,
                                generatedBy: String = "ShapeIcpProposal"
                              )(
                                implicit rand: scalismo.utils.Random
                              ) extends ProposalGenerator[ModelFittingParameters]
  with TransitionProbability[ModelFittingParameters] {

  private val referenceMesh = model.referenceMesh
  private val cashedPosterior: Memoize[ModelFittingParameters, LowRankGaussianProcess[_3D, EuclideanVector[_3D]]] = Memoize(icpPosterior, 20)


  private lazy val interpolatedModel = model.gp.interpolate(NearestNeighborInterpolator())
  val pcaLandmark = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "/Ref/landmarks/ref_landmarks.json")).get
  private val modelPoints = MyLMSampler3D(model.referenceMesh,Seq(pcaLandmark(0).point), 1000,1300).sample().map(_._1)
  private val modelIds: IndexedSeq[PointId] = modelPoints.map(p => model.referenceMesh.pointSet.findClosestPoint(p).id).toIndexedSeq
  private val targetPoints: IndexedSeq[Point[_3D]] = UniformMeshSampler3D(target, numOfSamplePoints).sample().map(_._1).toIndexedSeq

  override def propose(theta: ModelFittingParameters): ModelFittingParameters = {
    val posterior = cashedPosterior(theta)
    val proposed: Field[_3D, EuclideanVector[_3D]] = posterior.sample()

    def f(pt: Point[_3D]): Point[_3D] = pt + proposed(pt)

    val newCoefficients = model.coefficients(referenceMesh.transform(f))

    val currentShapeCoefficients = theta.shapeParameters.parameters
    val newShapeCoefficients = currentShapeCoefficients + (newCoefficients - currentShapeCoefficients) * stepLength

    theta.copy(
      shapeParameters = ShapeParameters(newShapeCoefficients),
      generatedBy = generatedBy
    )
  }


  private def randomStepLength(theta: ModelFittingParameters): Double = {
    1 - scala.util.Random.nextDouble() * 2.0
  }


  override def logTransitionProbability(from: ModelFittingParameters, to: ModelFittingParameters): Double = {
    val posterior = cashedPosterior(from)

    val compensatedTo = to.copy(shapeParameters = ShapeParameters(from.shapeParameters.parameters + (to.shapeParameters.parameters - from.shapeParameters.parameters) / stepLength))

    val pdf = posterior.logpdf(compensatedTo.shapeParameters.parameters)
    pdf
  }

  private def icpPosterior(theta: ModelFittingParameters): LowRankGaussianProcess[_3D, EuclideanVector[_3D]] = {
    def modelBasedClosestPointsEstimation(
                                           currentMesh: TriangleMesh[_3D],
                                           inversePoseTransform: RigidTransformation[_3D]
                                         ): IndexedSeq[(Point[_3D], EuclideanVector[_3D], MultivariateNormalDistribution)] = {

      //      val modelPoints = UniformMeshSampler3D(model.referenceMesh, numOfSamplePoints).sample().map(_._1)
      //      val modelIds: IndexedSeq[PointId] = modelPoints.map(p => model.referenceMesh.pointSet.findClosestPoint(p).id)
    //  val pointIds = target.pointSet.points.toIndexedSeq.map{
   //     pt => currentMesh.pointSet.findClosestPoint(currentMesh.operations.closestPointOnSurface(pt).point).id
   //   }
      val currentPoints = modelIds.map(id => (id, currentMesh.pointSet.point(id)))

      val noisyCorrespondence = currentPoints.map {
        case (id, pt) =>

          val closesttargetPoint = target.operations.closestPointOnSurface(pt).point
          val intersection = target.operations.getIntersectionPoints(pt,currentMesh.vertexNormals.atPoint(id) )

          val pointOnNormal =  if(intersection != Empty) {
            intersection.minBy { intersectionTargetPoint =>
              (intersectionTargetPoint - pt).norm2
            }
          } else Point(0,0,0)
          val targetMeshId = target.pointSet.findClosestPoint(pointOnNormal).id
          val normal = currentMesh.vertexNormals.atPoint(id)
          val normal2 = target.vertexNormals.atPoint(targetMeshId)
          val dot = (normal).dot(normal2)

          val flag =  if ((pointOnNormal - pt).norm2 < 1.0 && pointOnNormal!=Point(0,0,0) && ( dot <= 0))  true else false


          val targetPointId = target.pointSet.findClosestPoint(closesttargetPoint).id
          val isOnBoundary = currentMesh.operations.pointIsOnBoundary(id)
          val noiseDistribution = SurfaceNoiseHelpers.surfaceNormalDependantNoise(currentMesh.vertexNormals.atPoint(id), noiseAlongNormal, tangentialNoise)
          (id, pointOnNormal, noiseDistribution,   flag && isOnBoundary)
      }

      val correspondenceFiltered = if (boundaryAware) noisyCorrespondence.filter(!_._4) else noisyCorrespondence

      for ((pointId, targetPoint, uncertainty, _) <- correspondenceFiltered) yield {
        val referencePoint = model.referenceMesh.pointSet.point(pointId)
        (referencePoint, inversePoseTransform(targetPoint) - referencePoint, uncertainty)
      }
    }


    def targetBasedClosestPointsEstimation(
                                            currentMesh: TriangleMesh[_3D],
                                            inversePoseTransform: RigidTransformation[_3D]
                                          ): IndexedSeq[(Point[_3D], EuclideanVector[_3D], MultivariateNormalDistribution)] = {

      //      val targetPoints: Seq[Point[_3D]] = UniformMeshSampler3D(target, numOfSamplePoints).sample().map(_._1).toIndexedSeq

      val noisyCorrespondence = targetPoints.map { targetPoint =>
        val id = target.pointSet.findClosestPoint(targetPoint).id
        val intersection = currentMesh.operations.getIntersectionPoints(targetPoint,target.vertexNormals.atPoint(id) )

        val pointOnNormal =  if(intersection != Empty) {
          intersection.minBy { intersectionModelPoint =>
            (intersectionModelPoint - targetPoint).norm2
          }
        } else Point(0,0,0)

        val instanceId = currentMesh.pointSet.findClosestPoint(pointOnNormal).id
        val normal = currentMesh.vertexNormals.atPoint(instanceId)
        val normal2 = target.vertexNormals.atPoint(id)
        val dot = (normal).dot(normal2)

        val flag = if ((pointOnNormal - targetPoint).norm2 < 1.0 && pointOnNormal!=Point(0,0,0) && ( dot <= 0) ) true else false

        val isOnBoundary = target.operations.pointIsOnBoundary(id)
        val noiseDistribution = SurfaceNoiseHelpers.surfaceNormalDependantNoise(target.vertexNormals.atPoint(id), noiseAlongNormal, tangentialNoise)
        (instanceId, targetPoint, noiseDistribution,  flag && isOnBoundary)
      }

      val correspondenceFiltered = if (boundaryAware) noisyCorrespondence.filter(!_._4) else noisyCorrespondence

      for ((pointId, targetPoint, uncertainty, _) <- correspondenceFiltered) yield {
        val referencePoint = model.referenceMesh.pointSet.point(pointId)
        // (reference point, deformation vector in model space starting from reference, usually zero-mean observation uncertainty)
        (referencePoint, inversePoseTransform(targetPoint) - referencePoint, uncertainty)
      }
    }

    /**
      * Estimate where points should move to together with a surface normal dependant noise.
      *
      * @param theta Current fitting parameters
      * @return List of points, with associated deformation and normal dependant surface noise.
      */
    def uncertainDisplacementEstimation(theta: ModelFittingParameters)
    : IndexedSeq[(Point[_3D], EuclideanVector[_3D], MultivariateNormalDistribution)] = {
      val currentMesh = ModelFittingParameters.transformedMesh(model, theta)
      val inversePoseTransform = ModelFittingParameters.poseTransform(theta).inverse

      if (projectionDirection == TargetSampling) {
        targetBasedClosestPointsEstimation(currentMesh, inversePoseTransform)
      } else {
        modelBasedClosestPointsEstimation(currentMesh, inversePoseTransform)
      }
    }

    val uncertainDisplacements = uncertainDisplacementEstimation(theta)
    interpolatedModel.posterior(uncertainDisplacements)
  }

}
case class MyLMSampler3D(mesh: TriangleMesh[_3D], lmPoints: Seq[Point[_3D]], numOfNeightbours: Int, numberOfPoints: Int)(implicit rng: Random)
  extends Sampler[_3D] {

  override val volumeOfSampleRegion = mesh.area

  private val p: Double = 1.0 / mesh.area

  val samplePointsInit: Set[Point[_3D]] = lmPoints.flatMap(lm => (mesh.pointSet.findNClosestPoints(lm, numOfNeightbours).map(_.point))).toSet
  val samplePoints: IndexedSeq[(Point[_3D], Double)] = scala.util.Random.shuffle(samplePointsInit).take(numberOfPoints).map(f => (f, p)).toIndexedSeq



  override def sample() = samplePoints
}
