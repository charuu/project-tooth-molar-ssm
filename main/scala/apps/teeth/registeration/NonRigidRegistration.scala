package apps.teeth.registeration

import java.io.File

import apps.teeth.utilities.{LoadTestData, MyLMSampler3D, alignShapes}
import apps.teeth.utilities.MyLMSampler3D.ui
import apps.teeth.utilities.Paths.generalPath
import apps.util.AlignmentTransforms
import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common.{NearestNeighborInterpolator, PointId}
import scalismo.geometry.{EuclideanVector, Landmark, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.registration.GaussianProcessTransformationSpace
import scalismo.statisticalmodel.{MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless, ShowInScene, StatisticalMeshModelViewControls}
import scalismo.ui.model

object NonRigidRegistration {
  implicit val rng = scalismo.utils.Random(42)
  val ui = ScalismoUI()

  def main(args: Array[String]): Unit = {

    scalismo.initialize()

    val (model, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget("Partial")
    val alignedMeshLms = alignShapes(modelLms,targetMeshes,targetLms).align()

    register(model,modelLms,alignedMeshLms)

  }

  def register(model: StatisticalMeshModel, modelLms: Seq[Landmark[_3D]], alignedMeshLms: Seq[((TriangleMesh[_3D],String),(Seq[Landmark[_3D]],String))]) = {
    implicit val rng = scalismo.utils.Random(42)

    val varianceSeq = Seq(1.0,0.1, 0.05, 0.01)
    val varianceSeq2 = Seq(0.05)

    alignedMeshLms.take(1).map { targetMeshLms =>

      val target = targetMeshLms._1._1
      val targetName = targetMeshLms._1._2
      val targetLms = targetMeshLms._2._1

      val modelGroup = ui.createGroup("model")
      val showModel = ui.show(modelGroup, model, "model")

      val targetGroup = ui.createGroup("target")
      ui.show(targetGroup, target, "target")

      val TargetToModelFit = nonrigidICP( model,modelLms, target, targetLms, showModel, 30, varianceSeq, "TargetToModel")
      val pars = showModel.shapeModelTransformationView.shapeTransformationView.coefficients
      model.instance(pars)
      MeshIO.writeVTK(model.instance(pars), new File(generalPath + s"/Registered/mesh/TargetToModelFit.vtk"))

      val ModelToTargetFit = nonrigidICP(TargetToModelFit,modelLms, target, targetLms, showModel, 30, varianceSeq2, "ModelToTarget")
      val pars2 = showModel.shapeModelTransformationView.shapeTransformationView.coefficients
      MeshIO.writeVTK(model.instance(pars2), new File(generalPath + s"/Registered/mesh/ModelToTargetFit.vtk"))

      val ad = MeshMetrics.avgDistance(model.instance(pars2), target)
      println(s"avgDistance ${ad}")
    }

  }



  def nonrigidICP( movingMesh: StatisticalMeshModel, modelLm: Seq[Landmark[_3D]], targetMesh: TriangleMesh[_3D], targetLm:Seq[Landmark[_3D]], showModel:ShowInScene.ShowInSceneStatisticalMeshModel.View, numberOfIterations: Int, variance:Seq[Double], Mode:String): StatisticalMeshModel = {
    val commonLmNames = modelLm.map(_.id) intersect targetLm.map(_.id)
    val corr = commonLmNames.map(name => (movingMesh.referenceMesh.pointSet.findClosestPoint(modelLm.find(_.id == name).get.point).id, targetLm.find(_.id == name).get.point, true))
    val landmarkCorrFit = fitModel(corr, movingMesh, 0.01)

    showModel.shapeModelTransformationView.shapeTransformationView.coefficients = landmarkCorrFit
    ui.show(movingMesh.instance(landmarkCorrFit),"landmarkfit")

    variance.map{
      v =>
        ICP(corr, movingMesh ,modelLm , targetMesh, targetLm, showModel, 30 , v ,Mode)
    }.toIndexedSeq.last

  }

  def ICP(corr: Seq[(PointId, Point[_3D], Boolean)], movingMesh: StatisticalMeshModel,modelLms: Seq[Landmark[_3D]], targetMesh: TriangleMesh[_3D], targetLM: Seq[Landmark[_3D]], showModel: StatisticalMeshModelViewControls, i: Int, v: Double, Mode: String): StatisticalMeshModel = {
    if (i == 0) movingMesh

    else {
      val newcorr = if (Mode == "TargetToModel") {
        val sampler = UniformMeshSampler3D(targetMesh, 1000)
        val points: Seq[Point[_3D]] = sampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
        val ptIds = points.map(point => targetMesh.pointSet.findClosestPoint(point).id)

        attributeCorrespondences(movingMesh.mean,targetMesh, ptIds).filter(!_._3)

      } else {

        val lmNames = Seq("F1", "F2", "L1", "L2")
        val lmPoints = modelLms.filter(lm => lmNames.contains(lm.id)).map(_.point)
        val mySampler = MyLMSampler3D(movingMesh.mean, lmPoints, 5000, 1000)
        val ptIds = mySampler.sample().map { p => movingMesh.mean.pointSet.pointId(p._1).get }.toIndexedSeq

        attributeCorrespondences( movingMesh.mean,targetMesh, ptIds).filter(!_._3)
      }

      val pars = fitModel(newcorr, movingMesh, v)

      showModel.shapeModelTransformationView.shapeTransformationView.coefficients = pars
      println(s"iteration $i is ${i} ,variance ${v}")

      println(MeshMetrics.avgDistance(targetMesh, movingMesh.instance(pars)))
      ICP(newcorr, movingMesh, modelLms, targetMesh, targetLM, showModel, i - 1, v, Mode)
    }
  }


  def fitModel(correspondences: Seq[(PointId, Point[_3D],Boolean)],movingMesh:StatisticalMeshModel,v:Double): DenseVector[Double] = {
    val variance =  DenseMatrix.eye[Double](3) * v
    val littleNoise = MultivariateNormalDistribution(DenseVector.zeros[Double](3),variance )

    val regressionData = correspondences.map(correspondence =>
      (correspondence._1, correspondence._2, littleNoise)
    )
    val posterior = movingMesh.posterior(regressionData.toIndexedSeq)
    val pars = movingMesh.coefficients(posterior.mean)

    println(s"variance ${v}")

    pars
  }

  def attributeCorrespondences(movingMesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D], ptIds: Seq[(PointId)]): Seq[(PointId, Point[_3D], Boolean)] = {
    val findID = ptIds.map { id: PointId =>
      val closestPointOnTargetMesh = targetMesh.pointSet.point(id)
      val movingMeshId = movingMesh.pointSet.findClosestPoint(closestPointOnTargetMesh).id
      val normal = movingMesh.vertexNormals.atPoint(movingMeshId)
      val normal2 = targetMesh.vertexNormals.atPoint(id)
      val dot = normal2.dot(normal)
      //Math.round(dot)==1
      val getId = if (math.round(dot) >= 0) {
        true
      } else {
        false
      }
      (id, getId)
    }

    //if (p._2 == true)
    val x = for (p <- findID ) yield {

      val closestPointOnTargetMesh = targetMesh.pointSet.point(p._1)
      val movingmeshId = movingMesh.pointSet.findClosestPoint(closestPointOnTargetMesh).id
      val isOnBoundary = targetMesh.operations.pointIsOnBoundary(p._1)
      //  println((movingMesh.pointSet.findClosestPoint(closestPointOnMesh).point - closestPointOnMesh).norm,"distance between landmarks")

      (movingmeshId, closestPointOnTargetMesh, isOnBoundary)
    }
    println(s"count of points  : ${x.length}")
    x
  }

}
