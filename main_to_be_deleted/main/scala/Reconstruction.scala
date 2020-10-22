package main.scala
import java.io.File

import api.sampling.SurfaceNoiseHelpers
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.LoadTestData
import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common.{PointId, RealSpace}
import scalismo.geometry.{EuclideanVector, Landmark, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel, PDKernel}
import scalismo.mesh.{MeshMetrics, TriangleMesh}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.registration.LandmarkRegistration
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}
import scalismo.ui.api.ShowInScene.ShowInSceneStatisticalMeshModel.View

object Reconstruction {
  val ui = ScalismoUIHeadless()
  val path = "data/LowerMolar/"
  val (modl, modelLms, targetMesh, targetLms) = LoadTestData.modelAndTarget2()
  val targetGroup = ui.createGroup("Target Group")
  val SSMGroup = ui.createGroup("SSM Group")
  val middlept = LandmarkIO.readLandmarksJson[_3D](new File(generalPath + "/Registered/model/pcaModel_1000points_Landmarks_middle.json")).get
  val id3= modl.mean.pointSet.pointId(middlept.seq(0).point).get
  implicit val rng = scalismo.utils.Random(42)
 val movingmesh = modl
  val newModelGroup = ui.createGroup("New Result")
  val showModel = ui.show(newModelGroup,movingmesh,"Model")

  def main(args: Array[String]): Unit = {
   scalismo.initialize()

// Load cropped model
  val croppedMesh = targetMesh// MeshIO.readMesh(new File(generalPath + "/aligned/partial/mesh/" + "0012_36.stl")).get
  ui.show(targetGroup,croppedMesh,"cropped Mesh")

  val croppedLandmark = targetLms// LandmarkIO.readLandmarksJson[_3D](new File(generalPath + "/aligned/partial/landmarks/0012_36.json")).get
  val referenceLandmarkViews2 = croppedLandmark.map(lm => ui.show( targetGroup,lm, s"lm-${lm.id}"))
   // Load SSM
    val modelLandmark = modelLms//LandmarkIO.readLandmarksJson[_3D](new File(generalPath + "/Registered/model/AugCkModelLandmarks_500.json")).get
     val referenceLandmarkViews = modelLandmark.map(lm => ui.show( SSMGroup,lm, s"lm-${lm.id}"))

   nonRigidMCMC(modl.mean,croppedMesh)
 }

  def nonRigidMCMC(aligned: TriangleMesh[_3D],cropped: TriangleMesh[_3D]):TriangleMesh[_3D]= {
    implicit val rng = scalismo.utils.Random(42)
    val i =20
  //  val m = Seq(0.4,0.3,0.1,0.01,0.001,0.0001).map { v =>
    val m1 =  nonRigidICP(  aligned, cropped, i,0.4, showModel)
      val m2 =  nonRigidICP(  m1, cropped, i,0.3, showModel)
      val m3 =  nonRigidICP(  m2, cropped, i,0.1, showModel)
      val m4 =  nonRigidICP(  m3, cropped, i,0.01, showModel)
      val m5 =  nonRigidICP(  m4, cropped, i,0.001, showModel)
    val m6 =  nonRigidICP(  m5, cropped, i,0.0001, showModel)
    val m7 =  nonRigidICP(  m6, cropped, i,0.00001, showModel)
  //  }
    m7

  }

  def nonRigidICP(model: TriangleMesh[_3D], target: TriangleMesh[_3D], numberOfIterations: Int ,v:Double,showModel :View ) : TriangleMesh[_3D] = {

    if (numberOfIterations == 0) {
      val modelGroup = ui.createGroup("Result2")

      ui.show(modelGroup, model, "target Posterior")
      model}
    else
  {
        val targetSampler = UniformMeshSampler3D(target, numberOfPoints = 4000)
        val targetPoints: Seq[Point[_3D]] = targetSampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
        val targetPtIds = targetPoints.map(point => target.pointSet.findClosestPoint(point).id)

        val correspondences = attributeCorrespondencesTargetToModel(model, target, targetPtIds)
        val correspondenceFiltered = correspondences.filter(!_._3)

        val regressionTargetSampledData = correspondenceFiltered.map { correspondence =>
          val noiseDistribution = SurfaceNoiseHelpers.surfaceNormalDependantNoise(model.vertexNormals.atPoint(correspondence._1), v, v)

          (correspondence._1, correspondence._2, noiseDistribution)
        }

        val posteriorFromTarget = movingmesh.posterior(regressionTargetSampledData.toIndexedSeq)
        showModel.shapeModelTransformationView.shapeTransformationView.coefficients = movingmesh.coefficients(posteriorFromTarget.mean)
        println("applied posterior from target sampled data")
        val topToothPtIDs = posteriorFromTarget.mean.pointSet.pointIds.filter { id =>

          (posteriorFromTarget.mean.pointSet.point(id) - posteriorFromTarget.mean.pointSet.point(id3)).norm <= 5.5
        }

        val topToothPtIDsFromPosteriorModel = posteriorFromTarget.marginal(topToothPtIDs.toIndexedSeq)
        //  ui.show(topToothPtIDsFromPosteriorModel.mean,"top")

        val sampler = UniformMeshSampler3D(topToothPtIDsFromPosteriorModel.mean, numberOfPoints = 4000)
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

        showModel.shapeModelTransformationView.shapeTransformationView.coefficients = movingmesh.coefficients(posteriorFromModel.mean)
        println("applied posterior from target and model sampled data")

        println(MeshMetrics.avgDistance(target, posteriorFromModel.mean))

        println(s"iteration $numberOfIterations is ${numberOfIterations}, variance ${v}")
    //  val modelGroup = ui.createGroup("Result")

    //  ui.show(modelGroup, posteriorFromModel.mean, "target Posterior")
     // model.copy(posteriorFromModel.mean)
      nonRigidICP(posteriorFromModel.mean, target,  numberOfIterations -1 ,v,  showModel)


      }



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
    println(s"count of points  : ${count}")
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
    println(s"count of points  : ${count}")
    x
  }


}
