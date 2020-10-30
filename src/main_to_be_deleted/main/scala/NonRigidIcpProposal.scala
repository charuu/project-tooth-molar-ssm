package main.scala
import main.scala.ModelFittingParameters
import breeze.linalg.{DenseMatrix, DenseVector}
import main.scala.Reconstruction.{id3, movingmesh, ui}
import scalismo.common.PointId
import scalismo.geometry.{Point, _3D}
import scalismo.mesh.{MeshMetrics, TriangleMesh, TriangleMesh3D}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.sampling.{ProposalGenerator, TransitionProbability}
import scalismo.statisticalmodel.{MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.ShowInScene.ShowInSceneStatisticalMeshModel.View

case class NonRigidIcpProposal(
                                model: StatisticalMeshModel,
                                target: TriangleMesh3D,
                                variance: Double,
                                numOfSamplePoints: Int,
                                projectionDirection: String = "TargetSampling",
                                boundaryAware: Boolean = true,
                                generatedBy: String = "ShapeIcpProposal"
                              )(
                                implicit rand: scalismo.utils.Random
                              ) extends ProposalGenerator[ModelFittingParameters]
  with TransitionProbability[ModelFittingParameters] {
  val sampler = UniformMeshSampler3D(target, numberOfPoints = numOfSamplePoints )
  val points: Seq[Point[_3D]] = sampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
  val ptIds = points.map(point => target.pointSet.findClosestPoint(point).id)
  val finalFit = model.mean
  val corr2 = attributeCorrespondences(finalFit,target, ptIds)

  val perturbationDistr = new MultivariateNormalDistribution(
    DenseVector.zeros(model.rank),
    DenseMatrix.eye[Double](model.rank) * variance * variance
  )


  override def propose(sample: ModelFittingParameters): ModelFittingParameters = {
    val pars = nonrigidICP(model.mean, target, numOfSamplePoints, variance)
    val newParameters = sample.shapeParameters.copy(parameters = model.coefficients(pars.mean) + perturbationDistr.sample)
    sample.copy(sample.scalaParameter,sample.poseParameters,shapeParameters = newParameters ,generatedBy = s"ShapeUpdateProposal ($variance)")
  }

  def logTransitionProbability(from: ModelFittingParameters, to: ModelFittingParameters) = {
    val residual = to.shapeParameters.parameters - from.shapeParameters.parameters
    perturbationDistr.logpdf(residual)
  }

  def attributeCorrespondences(movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D],ptIds:Seq[(PointId)] ): Seq[(PointId, Point[_3D])]= {
    val findID = ptIds.map { id: PointId =>
      val targetPt = targetMesh.pointSet.point(id)
      val modelId = movingMesh.pointSet.findClosestPoint(targetPt).id
      val getId = if (math.round(movingMesh.vertexNormals.atPoint(modelId).dot(targetMesh.vertexNormals.atPoint(id))) == (1)) { true } else {  false  }
      (id, getId)
    }
    val x = for (p <- findID if (p._2 == true)) yield {
      val pt = targetMesh.pointSet.point(p._1)
      val closestPointOnMesh2 = movingMesh.pointSet.findClosestPoint(pt).id
      (closestPointOnMesh2,pt)
    }
    x
  }
  def fitModel(correspondences: Seq[(PointId, Point[_3D])],v:Double): TriangleMesh[_3D] = {
    val variance =  DenseMatrix.eye[Double](3) * v
    val littleNoise = MultivariateNormalDistribution(DenseVector.zeros[Double](3),variance )

    val regressionData = correspondences.map(correspondence =>
      (correspondence._1, correspondence._2, littleNoise)
    )
    val posterior = model.posterior(regressionData.toIndexedSeq)
    posterior.mean
  }
  /*def nonrigidICP(corr:Seq[(PointId,Point[_3D])],movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D], ptIds: Seq[PointId], numberOfIterations: Int,variance:Double): DenseVector[Double] = {
    if (numberOfIterations == 0) model.coefficients(movingMesh)
    else {
      val transformed = fitModel(corr, variance)
      nonrigidICP(attributeCorrespondences(movingMesh, targetMesh, ptIds),transformed, targetMesh,ptIds, numberOfIterations - 1,variance)
    }
  }*/

  def nonrigidICP(model: TriangleMesh[_3D], target: TriangleMesh[_3D], numberOfIterations: Int ,v:Double ) : StatisticalMeshModel = {


      val targetSampler = UniformMeshSampler3D(target, numberOfPoints = numOfSamplePoints)
      val targetPoints: Seq[Point[_3D]] = targetSampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
      val targetPtIds = targetPoints.map(point => target.pointSet.findClosestPoint(point).id)

      val correspondences = attributeCorrespondencesTargetToModel(model, target, targetPtIds)
      val correspondenceFiltered = correspondences.filter(!_._3)

      val regressionTargetSampledData = correspondenceFiltered.map { correspondence =>
        val noiseDistribution = SurfaceNoiseHelpers.surfaceNormalDependantNoise(model.vertexNormals.atPoint(correspondence._1), v, v)

        (correspondence._1, correspondence._2, noiseDistribution)
      }

      val posteriorFromTarget = movingmesh.posterior(regressionTargetSampledData.toIndexedSeq)
   //   showModel.shapeModelTransformationView.shapeTransformationView.coefficients = movingmesh.coefficients(posteriorFromTarget.mean)
   //   println("applied posterior from target sampled data")
      val topToothPtIDs = posteriorFromTarget.mean.pointSet.pointIds.filter { id =>
        (posteriorFromTarget.mean.pointSet.point(id) - posteriorFromTarget.mean.pointSet.point(id3)).norm <= 5.5
      }

      val topToothPtIDsFromPosteriorModel = posteriorFromTarget.marginal(topToothPtIDs.toIndexedSeq)
      //  ui.show(topToothPtIDsFromPosteriorModel.mean,"top")

      val sampler = UniformMeshSampler3D(topToothPtIDsFromPosteriorModel.mean, numberOfPoints = numOfSamplePoints)
      val modelPoints: Seq[Point[_3D]] = sampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
      val modelPtIds = modelPoints.map(point => posteriorFromTarget.mean.pointSet.findClosestPoint(point).id)
      var correspondences2 = attributeCorrespondencesModelToTarget(posteriorFromTarget.mean, target, modelPtIds)

      val correspondenceFiltered2 = correspondences2.filter(!_._3)
      correspondences2 = correspondences ++ correspondenceFiltered2
      val regressionData2 = correspondences2.map { correspondence =>
        val noiseDistribution = SurfaceNoiseHelpers.surfaceNormalDependantNoise(posteriorFromTarget.mean.vertexNormals.atPoint(correspondence._1), v, v)

        (correspondence._1, correspondence._2, noiseDistribution)
      }

      val posteriorFromModel = posteriorFromTarget.posterior(regressionData2.toIndexedSeq)

    //  showModel.shapeModelTransformationView.shapeTransformationView.coefficients = movingmesh.coefficients(posteriorFromModel.mean)
  //    println("applied posterior from target and model sampled data")

      println(MeshMetrics.avgDistance(target, posteriorFromModel.mean))

     // println(s"iteration $numberOfIterations is ${numberOfIterations}, variance ${v}")
      //  val modelGroup = ui.createGroup("Result")

      //  ui.show(modelGroup, posteriorFromModel.mean, "target Posterior")
      // model.copy(posteriorFromModel.mean)
     // nonRigidICP(posteriorFromModel.mean, target,  numberOfIterations -1 ,v,  showModel)

      posteriorFromModel




  }
  def attributeCorrespondencesModelToTarget(movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D],ptIds:Seq[(PointId)] ): Seq[(PointId, Point[_3D],Boolean)]= {
    val findID = ptIds.map { id: PointId =>
      val closestPointOnMesh = movingMesh.pointSet.point(id)
      val normal = movingMesh.vertexNormals.atPoint(id)
      val targetId = targetMesh.pointSet.findClosestPoint(closestPointOnMesh).id

      val normal2 = targetMesh.vertexNormals.atPoint(targetId)
      val dot = normal.dot(normal2)

      val getId = if (math.round(dot)>0) {true} else {false}
      (id, getId)
    }
    var count = 0
    val x = for (p <- findID if (p._2 == true)) yield {
      val isOnBoundary = movingMesh.operations.pointIsOnBoundary(p._1)
      val closestPointOnMesh = movingMesh.pointSet.point(p._1)
      val targetPoint = targetMesh.pointSet.findClosestPoint(closestPointOnMesh).point
      //  println((movingMesh.pointSet.findClosestPoint(closestPointOnMesh).point - closestPointOnMesh).norm,"distance between landmarks")
      count = count + 1
      (p._1,targetPoint,isOnBoundary)
    }
  //  println(s"count of points  : ${count}")
    x
  }

  def attributeCorrespondencesTargetToModel(movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D],ptIds:Seq[(PointId)] ): Seq[(PointId, Point[_3D],Boolean)]= {
    val findID = ptIds.map { id: PointId =>
      val closestPointOnTargetMesh = targetMesh.pointSet.point(id)
      val movingMeshId = movingMesh.pointSet.findClosestPoint(closestPointOnTargetMesh).id
      val normal = movingMesh.vertexNormals.atPoint(movingMeshId)
      val normal2 = targetMesh.vertexNormals.atPoint(id)
      val dot = normal.dot(normal2)
      //Math.round(dot)==1
      val getId = if (math.round(dot)>0) {true} else {false}
      (id, getId)
    }
    var count = 0
    val x = for (p <- findID if (p._2 == true)) yield {

      val closestPointOnTargetMesh = targetMesh.pointSet.point(p._1)
      val movingmeshId = movingMesh.pointSet.findClosestPoint(closestPointOnTargetMesh).id
      val isOnBoundary = movingMesh.operations.pointIsOnBoundary(movingmeshId)
      //  println((movingMesh.pointSet.findClosestPoint(closestPointOnMesh).point - closestPointOnMesh).norm,"distance between landmarks")
      count = count + 1
      (movingmeshId, closestPointOnTargetMesh,isOnBoundary)
    }
  //  println(s"count of points  : ${count}")
    x
  }
}