import java.io.File

import breeze.linalg.DenseMatrix
import changepointKernel.ui
import scalismo.common.{Field, NearestNeighborInterpolator, RealSpace}
import scalismo.geometry.{EuclideanVector, Landmark, Point, _3D}
import scalismo.io.StatisticalModelIO
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI

object changepointKernel {
  val ui = ScalismoUI()
  val lmList = new Array[Landmark[_3D]](400)
  val pcaModel = StatisticalModelIO.readStatisticalMeshModel(new java.io.File("Registered/model/pcaModel.h5")).get

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    implicit val rng = scalismo.utils.Random(42)
     val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File("Ref/model/Lower_500-components.h5")).get
    val gpSSM = model.gp.interpolate(NearestNeighborInterpolator())

    val zeroMean = Field(RealSpace[_3D], (pt:Point[_3D]) => EuclideanVector(0,0,0))
    val gk1 = DiagonalKernel(GaussianKernel[_3D](0.000001), 3)
    val gk11 = DiagonalKernel(GaussianKernel[_3D](0.0001), 3)
    val gk111 = DiagonalKernel(GaussianKernel[_3D](0.1), 3)
    val gk2 = DiagonalKernel(GaussianKernel[_3D](1.0), 3)

    val pt = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(20.6829776763916,27.43699836730957,0.1526535004377365)).point
    val pt2 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(17.485445022583008,25.554750442504883,2.171308994293213)).point
    val pt3 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(23.5267276763916,29.311439514160156,1.2742375135421753)).point
    val pt4 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(18.527441024780273,30.3780460357666,1.8765640258789062)).point
    val pt5 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(21.739791870117188,23.814760208129883,2.4957005977630615)).point
    val pt6 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(26.076618194580078,26.231367111206055,-1.2245844602584839)).point
    val pt7 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(21.734834671020508,32.82468032836914,-2.6155548095703125)).point
    val pt8 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(21.739791870117188,23.814760208129883,2.4957005977630615)).point
    val pt9 = pcaModel.referenceMesh.pointSet.findClosestPoint(Point(16.2012939453125,28.826005935668945,-2.062349319458008)).point

//16.2012939453125,28.826005935668945,-2.062349319458008
    ui.show(pcaModel.referenceMesh ,"pcaModel")

    val changePointKernel1 = ChangePointKernel(gk111, gk2,pt,1)
    val changePointKernel = ChangePointKernel(gk1, gk2,pt,1)
    val changePointKernel2 = ChangePointKernel(gk1, gk2,pt2,2)
    val changePointKernel3 = ChangePointKernel(gk1, gk2,pt3,3)
    val changePointKernel4 = ChangePointKernel(gk1, gk2,pt4,3)
    val changePointKernel5 = ChangePointKernel(gk1, gk2,pt5,3)
    val changePointKernel6 = ChangePointKernel(gk1, gk2,pt6,3)
    val changePointKernel7 = ChangePointKernel(gk1, gk2,pt7,3)
    val changePointKernel8 = ChangePointKernel(gk1, gk2,pt8,3)
    val changePointKernel9 = ChangePointKernel(gk1, gk2,pt9,3)
    val gpCP = GaussianProcess(zeroMean,  changePointKernel + changePointKernel + changePointKernel2 + changePointKernel3 + changePointKernel4 + changePointKernel5 + changePointKernel6 + changePointKernel7 + changePointKernel8  + changePointKernel9 )

    val sampleCP =  gpCP.sampleAtPoints(model.referenceMesh.pointSet)
    val sampleGroup = ui.createGroup("Sample")
    ui.show(sampleGroup, sampleCP, "ChangePointKernelGP_sample")
    val lowRankAugmentedGP = LowRankGaussianProcess.approximateGPCholesky(
      model.referenceMesh.pointSet,
      gpCP,
      relativeTolerance = 0.01,
      interpolator = NearestNeighborInterpolator(),
    )

    val augmentedSSM = StatisticalMeshModel(pcaModel.referenceMesh, lowRankAugmentedGP)
    ui.show(augmentedSSM,"SSM")
    StatisticalModelIO.writeStatisticalMeshModel(augmentedSSM,new File("/home/jain0000/IdeaProjects/icp-proposal/data/teeth/LowerMolar/Registered/model/AugmentedPCA.h5"))

  }
  case class ChangePointKernel(ker1 : MatrixValuedPDKernel[_3D], ker2: MatrixValuedPDKernel[_3D], center: Point[_3D], radius: Double) extends MatrixValuedPDKernel[_3D]{
    override val outputDim = 3
    def s(p: Point[_3D]): Double = {
      val distance = (center.toVector-p.toVector).norm
      if (distance > 2.0) {
       1
      }
       else (2)
    }
    def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
      val sx = s(x)
      val sy = s(y)
      val sx_id = pcaModel.referenceMesh.pointSet.pointId(x).get
      val sy_id = pcaModel.referenceMesh.pointSet.pointId(y).get
      pcaModel.gp.cov(sx_id,sy_id) + ker1(x,y) * (1-sx) * (1-sy)
    }
    override def domain: RealSpace[_3D] = RealSpace[_3D]
  }
}
