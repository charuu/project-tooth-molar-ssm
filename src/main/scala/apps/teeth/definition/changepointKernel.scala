package apps.teeth.definition

import java.io.File

import apps.teeth.utilities.Paths.generalPath
import breeze.linalg.DenseMatrix
import scalismo.common.{Field, RealSpace}
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.io.{LandmarkIO, StatisticalModelIO}
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI

object changepointKernel {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUI()
    implicit val rng = scalismo.utils.Random(42)

    val pcaLandmark = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "/Registered/model/pcaModel_1000points_Landmarks_middle.json")).get
    ui.show(pcaLandmark(0), "pcalm")

    val zeroMean = Field(RealSpace[_3D], (pt: Point[_3D]) => EuclideanVector(0, 0, 0))

    val pcaModel = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Registered/model/pcaModel_1000points_1405.h5")).get
    ui.show(ui.createGroup("PCA"), pcaModel, "pca")
    Seq(3.0).map { (s) =>
        val sigma = s.toString.replace(".","")
      Seq(400).map { b =>
        //    radius = r
        Seq(1000).map { p =>
          println("Creating changepoint kernel")
          val aug = ui.createGroup(s"aug$s$b$p")
          print("sigma: ", s, "radius", 5,"basis ", b)
          val gk0 = DiagonalKernel(GaussianKernel[_3D](1.0), 3) * 0.1
          val gk1 = DiagonalKernel(GaussianKernel[_3D](s), 3) * 0.1
          val ck = CreateChangePointKernel(pcaModel, gk0 , pcaLandmark(0).point, 5)
          val gpCP = GaussianProcess(zeroMean, ck)
          val sampler = UniformMeshSampler3D(pcaModel.mean, numberOfPoints = p)

          println("Low rank GP approximation")
          val lowRankAugmentedGP = LowRankGaussianProcess.approximateGPNystrom(gpCP, sampler, numBasisFunctions = b)
          println("Creation of SSM")
          val pp= 1000
          val augmentedSSM = StatisticalMeshModel.augmentModel(pcaModel, lowRankAugmentedGP)
          StatisticalModelIO.writeStatisticalMeshModel(augmentedSSM, new File(generalPath + s"/Registered/model/AugmentedPCA_1405B.h5"))
          ui.show(aug, augmentedSSM, "SSM")

        }
      }
    }
  }
  case class CreateChangePointKernel(pca : StatisticalMeshModel, ker: MatrixValuedPDKernel[_3D], center: Point[_3D], radius: Double) extends MatrixValuedPDKernel[_3D]{
    override val outputDim = 3
    val ref = pca.mean
    def s(p: Point[_3D]): Double = {
      val d = Math.max(0, Math.sqrt(Math.pow(p.x-center.x, 2) + Math.pow(p.y-center.y, 2) + Math.pow(p.z-center.z, 2))-radius)
      1.0/(1+math.exp(d))
    }
    def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
      val sx = s(x)
      val sy = s(y)
      val sx_id = ref.pointSet.findClosestPoint(x).id
      val sy_id = ref.pointSet.findClosestPoint(y).id
      val pcaGP = pca.gp.cov(sx_id,sy_id)
      (ker(x,y)) * sx * sy //+ pcaGP * (1-sx) * (1-sy)
//      ker(x,y) * sx * sy // Activate ker in area only - no variance everywhere else
    }
    override def domain: RealSpace[_3D] = RealSpace[_3D]
  }
}
