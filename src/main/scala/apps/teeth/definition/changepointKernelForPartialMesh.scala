package apps.teeth.definition

import java.io.File

import apps.teeth.utilities.Paths.generalPath
import breeze.linalg.DenseMatrix
import scalismo.common.{Domain, Field, NearestNeighborInterpolator, RealSpace}
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.io.{LandmarkIO, StatisticalModelIO}
import scalismo.kernels.{DiagonalKernel, DiscreteMatrixValuedPDKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, UniformMeshSampler3D}
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}

import scala.collection.immutable.Stream.Empty

object changepointKernelForPartialMesh {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUIHeadless()
    implicit val rng = scalismo.utils.Random(42)
    val pcaLandmark = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "/Ref/landmarks/ref_landmarks.json")).get

    ui.show(pcaLandmark(4), "pcalm")

    val pcaModel = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Registered/model/pcaModel.h5")).get
    val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Ref/model/tooth_gp_model_A-components.h5")).get

    val zeroMean = Field(RealSpace[_3D], (pt: Point[_3D]) =>EuclideanVector(0, 0, 0))
    ui.show(ui.createGroup("PCA"), pcaModel, "pca")
    val varianceSeq = Seq((3.0, 2.0, 0.0))
    //Seq((3.0, 0.0, 0.0),(2.0, 0.0, 0.0),(4.0, 2.0, 0.0),(3.0, 1.0, 0.0),(2.0, 1.0, 0.0),(1.0, 0.0, 0.0))
    val scaleSeq = Seq((0.04,0.03, 0.0))
    //Seq((0.01,0.0, 0.0),(0.01,0.0, 0.0),(0.01,0.05, 0.0),(0.01,0.05, 0.0),(0.01,0.05, 0.0),(0.1,0.0, 0.0))

    println("Creating changepoint kernel")
    val aug = ui.createGroup(s"aug1_1000")
    varianceSeq.zip(scaleSeq).seq.foreach { x =>
      val aug = ui.createGroup(s"aug1_1000${x._1},${x._2}")

      val gk0 = DiagonalKernel(GaussianKernel[_3D](x._1._1),3) * x._2._1
      val gk3 = if(x._1._2!=0) {
        gk0 + DiagonalKernel(GaussianKernel[_3D](x._1._2), 3) * x._2._2
      } else gk0

      val gk1 = if(x._1._3!=0) {
        gk0 + gk3 + DiagonalKernel(GaussianKernel[_3D](x._1._3), 3) * x._2._3
      } else gk0 + gk3

      val gk2 = pcaModel.gp.cov
      //val gk4 = DiagonalKernel(GaussianKernel[_3D](10),3)* 0.01
      val ck = CreateChangePointKernel(pcaModel, model, gk1 , gk2, pcaLandmark(4).point, 5)

      val gpCP = GaussianProcess(zeroMean, ck)

      //   val sampler = FixedPointsUniformMeshSampler3D(model.referenceMesh, numberOfPoints = 30000)

      println("Low rank GP approximation")
      val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
      val lowRankAugmentedGP = LowRankGaussianProcess.approximateGPCholesky(model.referenceMesh.pointSet,gpCP,0.01,interpolator )
      println("Creation of SSM")
      val augmentedSSM = StatisticalMeshModel(model.referenceMesh, lowRankAugmentedGP)
      StatisticalModelIO.writeStatisticalMeshModel(augmentedSSM, new File(generalPath + s"/Registered/model/AugmentedPCAD.h5"))

      ui.show(aug, augmentedSSM, "SSM")
    }

  }
  case class CreateChangePointKernel(pca : StatisticalMeshModel,model : StatisticalMeshModel, ker: MatrixValuedPDKernel[_3D],ker2: DiscreteMatrixValuedPDKernel[_3D], center: Point[_3D], radius: Double) extends MatrixValuedPDKernel[_3D]{
    override val outputDim = 3
    //  val ref = model.mean
    def s(p: Point[_3D]): Double = {
      val d = Math.sqrt(Math.pow(p.x-center.x, 2) + Math.pow(p.y-center.y, 2) + Math.pow(p.z-center.z, 2))-radius
      1.0/(1+math.exp(d))
    }
    def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
      val sx = s(x)
      val sy = s(y)
      val sx_id = model.mean.pointSet.findClosestPoint(x).id
      val sy_id = model.mean.pointSet.findClosestPoint(y).id
      //  val pcaGP = pca.gp.cov(sx_id,sy_id)
      val k = DiagonalKernel(GaussianKernel[_3D](10),3) *0.001
      //println("***",pcaGP,"***")
      ker(x,y) * (sx) * (sy)  +   ker2(sx_id,sy_id)  +  k(x,y)// model.gp.cov(sx_id,sy_id) *0.005
      // pca.gp.cov(sx_id,sy_id)

    }
    override def domain: RealSpace[_3D] = RealSpace[_3D]
  }



}
