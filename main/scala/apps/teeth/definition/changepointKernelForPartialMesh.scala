package apps.teeth.definition

import java.io.File

import apps.teeth.utilities.Paths.generalPath
import breeze.linalg.DenseMatrix
import scalismo.common.{Field, RealSpace}
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.io.{LandmarkIO, StatisticalModelIO}
import scalismo.kernels.{DiagonalKernel, DiscreteMatrixValuedPDKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, UniformMeshSampler3D}
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}

object changepointKernelForPartialMesh {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUIHeadless()
    implicit val rng = scalismo.utils.Random(42)
    val pcaLandmark = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "/Ref/landmarks/ref_landmarks_7.json")).get

    ui.show(pcaLandmark(4), "pcalm")

    val pcaModel = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Registered/model/pcaModel_1000points_0406.h5")).get
    val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Ref/model/tooth_gp_model_50000-components.h5")).get
     val zeroMean = Field(RealSpace[_3D], (pt: Point[_3D]) =>EuclideanVector(0, 0, 0))
    ui.show(ui.createGroup("PCA"), pcaModel, "pca")

    println("Creating changepoint kernel")
    val aug = ui.createGroup(s"aug1_1000")
    val gk0 = DiagonalKernel(GaussianKernel[_3D](1.0), 3) * 0.05
    val gk4 = DiagonalKernel(GaussianKernel[_3D](0.5), 3) * 0.01
    val gk3 = DiagonalKernel(GaussianKernel[_3D](2.0), 3) * 0.01
    val gk1 =  DiagonalKernel(GaussianKernel[_3D](3.0), 3) * 0.1
    val gk2 =  pcaModel.gp.cov
    val ck = CreateChangePointKernel(pcaModel, model, gk0 + gk1 + gk3 + gk4,gk2, pcaLandmark(4).point, 5)

    val gpCP = GaussianProcess(zeroMean, ck)

    val sampler = FixedPointsUniformMeshSampler3D(model.mean, numberOfPoints = 1000)

    println("Low rank GP approximation")
    val lowRankAugmentedGP = LowRankGaussianProcess.approximateGPNystrom(gpCP, sampler, numBasisFunctions = 300)
    println("Creation of SSM")
    val augmentedSSM = StatisticalMeshModel(model.referenceMesh, lowRankAugmentedGP)
    StatisticalModelIO.writeStatisticalMeshModel(augmentedSSM, new File(generalPath + s"/Registered/model/AugmentedPCA_0906.h5"))
    ui.show(aug, augmentedSSM, "SSM")


  }
  case class CreateChangePointKernel(pca : StatisticalMeshModel,model : StatisticalMeshModel, ker: MatrixValuedPDKernel[_3D],ker2: DiscreteMatrixValuedPDKernel[_3D], center: Point[_3D], radius: Double) extends MatrixValuedPDKernel[_3D]{
    override val outputDim = 3
  //  val ref = model.mean
    def s(p: Point[_3D]): Double = {

      val d = Math.sqrt(Math.pow(p.x-center.x, 2) + Math.pow(p.y-center.y, 2) + Math.pow(p.z-center.z, 2))-radius

     // println("difference :",d, "return value :" ,  1.0/(1+math.exp(d)) )
      1.0/(1+math.exp(d))
    //   d
    }
    def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
      val sx = s(x)
      val sy = s(y)
      val sx_id = model.mean.pointSet.findClosestPoint(x).id
      val sy_id = model.mean.pointSet.findClosestPoint(y).id
      //val pcaGP = pca.gp.cov(sx_id,sy_id)
      //println("***",pcaGP,"***")
       ker(x,y) * (sx) * (sy) +   ker2(sx_id,sy_id)
     // pca.gp.cov(sx_id,sy_id)

    }
    override def domain: RealSpace[_3D] = RealSpace[_3D]
  }



}
