package main.scala
import java.io.File

import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common._
import scalismo.geometry._
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.statisticalmodel.MultivariateNormalDistribution
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless, ShowInScene}
import main.scala.Paths.generalPath

object Registration {
  scalismo.initialize()
  val ui = ScalismoUIHeadless()
  val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath +"/Ref/model/tooth_gp_model_1000-components.h5")).get
 // val showModel = ui.show(ui.createGroup("model"), model, "model")
  implicit val rng = scalismo.utils.Random(42)
  val modellandmark = LandmarkIO.readLandmarksJson[_3D](new java.io.File(generalPath + "/Ref/landmarks/ref_landmarks.json")).get

  def main(args: Array[String]): Unit = {
    val varianceSeq = Seq(1,0.5,0.1)
    // aligning
    val sampleLmFileList = new java.io.File(generalPath +"/aligned/landmarks").listFiles.sortBy(f => f.getName())
    val sampleMeshesFileList = new java.io.File(generalPath+"/aligned/meshes").listFiles.sortBy(f => f.getName())

    val sampler = UniformMeshSampler3D(model.referenceMesh, numberOfPoints = 500)
    val points: Seq[Point[_3D]] = sampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
    val ptIds = points.map(point => model.referenceMesh.pointSet.findClosestPoint(point).id)
    val finalFit = model.mean

    sampleMeshesFileList.par.foreach{f =>

      val target = MeshIO.readMesh(f).get
      val targetLMfile = sampleLmFileList.find(l => l.getName == f.getName.replace(".stl", ".json")).get
      val targetLM = LandmarkIO.readLandmarksJson[_3D](targetLMfile).get
      val modelGroup = ui.createGroup("model")
      val showModel = ui.show(modelGroup, model, "model")

      varianceSeq.foreach { v =>
        println(s"${v} variance on ${f.getName}")
        ui.show(target, s"${f.getName}")
        val corr = (modellandmark zip targetLM).map(f => (model.referenceMesh.pointSet.findClosestPoint(f._1.point).id, f._2.point))
        println("Non rigid icp for lm corr")
        val finalFit2 = nonrigidICP(corr, finalFit, target, ptIds,showModel, 1, v)
        println("Non rigid icp for randomly sampled points")
        val corr2 = attributeCorrespondences(finalFit2, target, ptIds)
        val finalFit3 = nonrigidICP(corr2, finalFit2, target, ptIds,showModel, 100, v)

        val ad = MeshMetrics.avgDistance(finalFit3, target)
        val hd = MeshMetrics.hausdorffDistance(finalFit3, target)

        println(s"hausdorffDistance ${hd}")
        println(s"avgDistance ${ad}")

        ui.show(finalFit3, "final fit3")
        MeshIO.writeSTL(finalFit3, new File(generalPath +s"/Registered/mesh/${f.getName}"))
      }
    }

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

  def attributeCorrespondences(movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D],ptIds:Seq[(PointId)] ): Seq[(PointId, Point[_3D])]= {
    val findID = ptIds.map { id: PointId =>
      val pt = movingMesh.pointSet.point(id)
      val alId = targetMesh.pointSet.findClosestPoint(pt).id
      val getId = if (math.round(movingMesh.vertexNormals.atPoint(id).dot(targetMesh.vertexNormals.atPoint(alId))) == (1)) { true } else {  false  }
      (id, getId)
    }
    val x = for (p <- findID if (p._2 == true)) yield {
      val pt = movingMesh.pointSet.point(p._1)
      val closestPointOnMesh2 = targetMesh.pointSet.findClosestPoint(pt).point
      (p._1, targetMesh.pointSet.findClosestPoint(closestPointOnMesh2).point)
    }
    x
  }

  def nonrigidICP(corr:Seq[(PointId,Point[_3D])],movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D], ptIds: Seq[PointId],showModel:ShowInScene.ShowInSceneStatisticalMeshModel.View, numberOfIterations: Int,variance:Double): TriangleMesh[_3D] = {
    if (numberOfIterations == 0) movingMesh
    else {
      val transformed = fitModel(corr, variance)
      val pars = model.coefficients(transformed)
      showModel.shapeModelTransformationView.shapeTransformationView.coefficients = pars
      println(s"iteration $numberOfIterations is ${numberOfIterations}")
      nonrigidICP(attributeCorrespondences(movingMesh, targetMesh, ptIds),transformed, targetMesh,ptIds,showModel, numberOfIterations - 1,variance)
    }
  }
}
