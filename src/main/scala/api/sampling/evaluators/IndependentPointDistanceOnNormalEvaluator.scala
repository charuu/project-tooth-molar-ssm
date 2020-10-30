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

package api.sampling.evaluators

import java.io.File

import api.sampling.ModelFittingParameters
import api.sampling.proposals.MyLMSampler3D
import apps.teeth.utilities.Paths.generalPath
import breeze.stats.distributions.ContinuousDistr
import scalismo.common.PointId
import scalismo.geometry.{Point, _3D}
import scalismo.io.LandmarkIO
import scalismo.mesh.TriangleMesh3D
import scalismo.numerics.UniformMeshSampler3D
import scalismo.sampling.DistributionEvaluator
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.utils.Random.implicits._

import scala.collection.immutable.Stream.Empty

case class IndependentPointDistanceOnNormalEvaluator(model: StatisticalMeshModel,
                                             targetMesh: TriangleMesh3D,
                                             likelihoodModel: ContinuousDistr[Double],
                                             evaluationMode: EvaluationMode,
                                             numberOfPointsForComparison: Int)
  extends DistributionEvaluator[ModelFittingParameters] with EvaluationCaching {
  val pcaLandmark = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "/Ref/landmarks/ref_landmarks.json")).get
  def getRandomPointsOnTarget: IndexedSeq[Point[_3D]] = {
    if (numberOfPointsForComparison >= targetMesh.pointSet.numberOfPoints) {
      targetMesh.pointSet.points.toIndexedSeq
    }
    else {
      UniformMeshSampler3D(targetMesh, numberOfPointsForComparison).sample().map(_._1)
    }
  }

  def getRandomPointIdsOnModel: IndexedSeq[PointId] = {
    if (numberOfPointsForComparison >= model.referenceMesh.pointSet.numberOfPoints) {
      model.referenceMesh.pointSet.pointIds.toIndexedSeq
    }
    else {

      val modelPoints = MyLMSampler3D(model.referenceMesh,Seq(pcaLandmark(0).point), 1000,500).sample().map(_._1)
      val modelIds: IndexedSeq[PointId] = modelPoints.map(p => model.referenceMesh.pointSet.findClosestPoint(p).id).toIndexedSeq
      modelIds
    //  UniformMeshSampler3D(model.referenceMesh, numberOfPointsForComparison).sample().map(_._1)
  //      .map(p => model.referenceMesh.pointSet.findClosestPoint(p).id)
    }
  }

  // Make sure not to oversample if the numberOfPointsForComparison is set higher than the points in the target or the model
  private val randomPointsOnTarget: IndexedSeq[Point[_3D]] = getRandomPointsOnTarget
  private val randomPointIdsOnModel: IndexedSeq[PointId] = getRandomPointIdsOnModel

  def distModelToTarget(modelSample: TriangleMesh3D): Double = {
    val pointsOnSample = randomPointIdsOnModel.map(modelSample.pointSet.point)
    val dists = for (pt <- pointsOnSample) yield {
      val modelPointId = modelSample.pointSet.findClosestPoint(pt).id
      val intersection = targetMesh.operations.getIntersectionPoints(pt,modelSample.vertexNormals.atPoint(modelPointId) )

      val pointOnNormal =if(intersection!=Empty) {intersection.minBy { intersectionTargetPoint =>
          (intersectionTargetPoint - pt).norm2
        }}else Point(0,0,0)
      val isOnBoundary = modelSample.operations.pointIsOnBoundary(modelPointId)
      //targetMesh.operations.closestPointOnSurface(pt).point
      val flag =  if ((pointOnNormal - pt).norm2 < 1.0 && pointOnNormal!=Point(0,0,0)) true else false
      (pointOnNormal,pt,flag&&isOnBoundary)
    }
    dists.filter(_._3==true).map{ tuple=>likelihoodModel.logPdf((tuple._1 - tuple._2).norm) }.sum


  }


  def distTargetToModel(modelSample: TriangleMesh3D): Double = {
    val dists = for (pt <- randomPointsOnTarget) yield {
      val targetPointId = targetMesh.pointSet.findClosestPoint(pt).id

        val intersection = modelSample.operations.getIntersectionPoints(pt, targetMesh.vertexNormals.atPoint(targetPointId))
        val pointOnNormal = if(intersection!=Empty) { intersection.minBy { intersectionModelPoint =>
            (intersectionModelPoint - pt).norm2
          }}else Point(0,0,0)
// likelihoodModel.logPdf((modelSample.operations.closestPointOnSurface(pt).point - pt).norm)
      val isOnBoundary = targetMesh.operations.pointIsOnBoundary(targetPointId)
      val flag =  if ((pointOnNormal - pt).norm2 < 1.0 && pointOnNormal!=Point(0,0,0)) true else false
      (pointOnNormal,pt,flag && isOnBoundary)
    }
    // likelihoodModel.logPdf((pointOnNormal - pt).norm)
   // dists.sum
    dists.filter(_._3==true).map{ tuple=>likelihoodModel.logPdf((tuple._1 - tuple._2).norm) }.sum
  }


  def computeLogValue(sample: ModelFittingParameters): Double = {

    val currentSample = ModelFittingParameters.transformedMesh(model, sample)
    val dist = evaluationMode match {
      case ModelToTargetEvaluation => distModelToTarget(currentSample)
      case TargetToModelEvaluation => distTargetToModel(currentSample)
      case SymmetricEvaluation => 0.5 * distModelToTarget(currentSample) + 0.5 * distTargetToModel(currentSample)
    }
    dist
  }
}



