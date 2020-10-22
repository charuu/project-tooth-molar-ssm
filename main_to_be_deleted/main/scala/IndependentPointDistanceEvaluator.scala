

package main.scala

import main.scala.ModelFittingParameters
import breeze.stats.distributions.ContinuousDistr
import scalismo.common.PointId
import scalismo.geometry.{Point, _3D}
import scalismo.mesh.TriangleMesh3D
import scalismo.numerics.UniformMeshSampler3D
import scalismo.sampling.DistributionEvaluator
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.utils.Random.implicits._

case class IndependentPointDistanceEvaluator(model: StatisticalMeshModel, targetMesh: TriangleMesh3D,
                                             likelihoodModel: ContinuousDistr[Double],
                                             evaluationMode: String,
                                             numberOfPointsForComparison: Int)
  extends DistributionEvaluator[ModelFittingParameters]  {

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
      UniformMeshSampler3D(model.referenceMesh, numberOfPointsForComparison).sample().map(_._1)
        .map(p => model.referenceMesh.pointSet.findClosestPoint(p).id)
    }
  }

  // Make sure not to oversample if the numberOfPointsForComparison is set higher than the points in the target or the model
  private val randomPointsOnTarget: IndexedSeq[Point[_3D]] = getRandomPointsOnTarget
  private val randomPointIdsOnModel: IndexedSeq[PointId] = getRandomPointIdsOnModel

  def distModelToTarget(modelSample: TriangleMesh3D): Double = {
    val pointsOnSample = randomPointIdsOnModel.map(modelSample.pointSet.point)
    val dists = for (pt <- pointsOnSample) yield {
      likelihoodModel.logPdf((targetMesh.operations.closestPointOnSurface(pt).point - pt).norm)
    }
    dists.sum
  }


  def distTargetToModel(modelSample: TriangleMesh3D): Double = {
    val dists = for (pt <- randomPointsOnTarget) yield {
      likelihoodModel.logPdf((modelSample.operations.closestPointOnSurface(pt).point - pt).norm)
    }
    dists.sum
  }


  def computeLogValue(sample: ModelFittingParameters): Double = {

    val currentSample = ModelFittingParameters.transformedMesh(model, sample)
    val dist = evaluationMode match {
      case "ModelToTargetEvaluation" => distModelToTarget(currentSample)
      case "TargetToModelEvaluation" => distTargetToModel(currentSample)
      case "SymmetricEvaluation" => 0.5 * distModelToTarget(currentSample) + 0.5 * distTargetToModel(currentSample)
    }
    dist
  }

  override def logValue(sample: ModelFittingParameters): Double = 0.0
}

