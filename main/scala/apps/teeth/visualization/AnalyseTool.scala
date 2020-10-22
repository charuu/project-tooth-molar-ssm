package apps.teeth.visualization

import apps.teeth.utilities.LoadTestData
import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common.{DiscreteField, NearestNeighborInterpolator, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector, EuclideanVector3D, Point, _3D}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.statisticalmodel.{LowRankGaussianProcess, MultivariateNormalDistribution}
import scalismo.ui.api.ScalismoUI

object AnalyseTool {
  def main(args: Array[String]):Unit = {
      scalismo.initialize()
    implicit val rng = scalismo.utils.Random(42)
     val ui = ScalismoUI(s"ScalismoUI")
    val targetGroup = ui.createGroup("Target")
    val targetinitGroup = ui.createGroup("Targetinit")
    val modelGroup = ui.createGroup("modelGroup")
    val (model, modelLms, targetMesh, targetLms) = LoadTestData.modelAndTarget("Partial")

    //check variance using landmark points
  //  ui.show(modelGroup,model.referenceMesh, "referenceMesh")
    ui.show(modelGroup,model, "mean")
 //   ui.show(targetinitGroup,targetInit, "targetInit")

    val referenceLandmarkViews = targetLms.map { lm =>
      ui.show(targetGroup, lm._1, "target transform")
    }
 //    targetLmInit.map { lm =>
 //     ui.show(targetinitGroup, lm, "target init ")
  //  }
    ui.show(targetGroup, targetMesh.seq(0)._1, "targetMesh")


      ui.show( modelGroup,modelLms, "model")


    val referenceLandmarkViews2 = modelLms.map { lm =>
      val lmId = model.referenceMesh.pointSet.findClosestPoint(lm.point).id
      val lmReference = model.referenceMesh.pointSet.point(lmId)
      val lmTipMean = model.mean.pointSet.point(lmId)
      val lmTipDomain = UnstructuredPointsDomain(IndexedSeq(lmReference))
      val lmTipDeformationAsSeq = IndexedSeq((lmTipMean - lmReference) )

      println((lmTipMean - lmReference).norm,"distance")
      //print distance of model lm to target lm
      // (1.3533618502053413,distance)
      // (1.9318538398207372,distance)
      // (1.2378333011496694,distance)
      // (1.418312155922957,distance)

      //   val lmTipDeformationAsSeq = IndexedSeq(targetMesh.vertexNormals.atPoint(lmId)*3)
      val lmDeformationField = DiscreteField[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]](lmTipDomain, lmTipDeformationAsSeq)

   //   val observationGroup = ui.createGroup("observation")
 //     ui.show(observationGroup, lmDeformationField, "lmTip")

    }


    val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
    val contGP = model.gp.interpolate(interpolator)
    val gp : LowRankGaussianProcess[_3D, EuclideanVector[_3D]] = model.gp.interpolate(NearestNeighborInterpolator())

    val sampler = UniformMeshSampler3D(targetMesh.seq(0)._1, numberOfPoints = 10)
    val points: Seq[Point[_3D]] = sampler.sample.map(pointWithProbability => pointWithProbability._1) // we only want the points
    val ptIds = points.map(point => targetMesh.seq(0)._1.pointSet.findClosestPoint(point).id)

    val center: EuclideanVector[_3D] = targetMesh.seq(0)._1.pointSet.points.map(_.toVector).reduce(_ + _) * 1.0 / targetMesh.seq(0)._1.pointSet.numberOfPoints.toDouble
    var sum :Double= 0.0
    ptIds.map { id =>
      val noise = MultivariateNormalDistribution(DenseVector.zeros[Double](3), DenseMatrix.eye[Double](3))

      val lmTipMean = targetMesh.seq(0)._1.pointSet.point(id)

      val lmReference = model.mean.pointSet.findClosestPoint(lmTipMean).point

      val lmTipDeformationAsSeq = IndexedSeq(((lmTipMean ) - (lmReference) ) )
       println(((lmTipMean ) - (lmReference) ).norm)
      val completePointDomain =
        UnstructuredPointsDomain(IndexedSeq(lmReference))
      println("Is it on boundary of fixed mesh ? , ",completePointDomain.isDefinedAt(lmTipMean))
      val lmDeformationField = DiscreteField[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]](completePointDomain, lmTipDeformationAsSeq)

      val observationGroup = ui.createGroup("observation")
      ui.show(observationGroup, lmDeformationField, "lmTip")
      sum = sum + (lmTipMean - lmReference).norm



    }
    println("average distance ",sum/2000)


  }
}
