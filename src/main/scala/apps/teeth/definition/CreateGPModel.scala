package apps.teeth.definition

import java.io.File

import apps.teeth.utilities.Paths.generalPath
import breeze.linalg.DenseMatrix
import breeze.linalg.svd.SVD
import scalismo.common.{Domain, NearestNeighborInterpolator, RealSpace, VectorField}
import scalismo.geometry.{EuclideanVector, Point, SquareMatrix, _3D}
import scalismo.io.{MeshIO, StatismoIO}
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.mesh.TriangleMesh3D
import scalismo.numerics.UniformMeshSampler3D
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}
import scalismo.utils.Random

object CreateGPModel {
  implicit val random: Random = Random(1024)

  def approxTotalVariance(gp: GaussianProcess[_3D, EuclideanVector[_3D]], evaluationDomain: TriangleMesh3D): Double = {
    val sampler = UniformMeshSampler3D(evaluationDomain, numberOfPoints = 65000)
    val points = sampler.sample.unzip._1
    val covValues = for (pt <- points) yield {
      gp.cov(pt, pt)
    }
    val variance = covValues.map(cm => cm(0, 0) + cm(1, 1) + cm(2, 2)).sum /// points.size
    variance
  }

  def getAxisOfMainVariance(mesh: TriangleMesh3D): DenseMatrix[Double] = {
    val N = 1.0 / mesh.pointSet.numberOfPoints
    val centerOfMass = (mesh.pointSet.points.foldLeft[EuclideanVector[_3D]](EuclideanVector(0f, 0f, 0f))((acc, e) => acc + e.toVector) * N).toPoint
    val cov = mesh.pointSet.points.foldLeft[SquareMatrix[_3D]](SquareMatrix.zeros)((acc, e) => acc + (e - centerOfMass).outer(e - centerOfMass)) * N
    val SVD(u, _, _) = breeze.linalg.svd(cov.toBreezeMatrix)
    u
  }

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val ui = ScalismoUIHeadless()

    Seq(108642).foreach { i =>
      val referenceMesh = MeshIO.readMesh(new File(generalPath, "Ref/meshes/lowermolar_LowerJaw_full.stl")).get

      val outputModelFile = new File(generalPath, s"Ref/model/tooth_gp_model_C-components.h5")

      println("Num of points in ref: " + referenceMesh.pointSet.numberOfPoints)

      val zeroMean = VectorField(RealSpace[_3D], (_: Point[_3D]) => EuclideanVector.zeros[_3D])

      val cov: MatrixValuedPDKernel[_3D] = new MatrixValuedPDKernel[_3D]() {
        private val kernels =  DiagonalKernel[_3D](GaussianKernel(10), 3) * 1.0 + DiagonalKernel[_3D](GaussianKernel(6), 3) * 0.1 + DiagonalKernel[_3D](GaussianKernel(4), 3) * 0.05

        override protected def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
          kernels(x, y)
        }

        override def outputDim = 3

        override def domain: Domain[_3D] = RealSpace[_3D]
      }

      val gp = GaussianProcess[_3D, EuclideanVector[_3D]](zeroMean, cov)

      val totalVariance = approxTotalVariance(gp, referenceMesh)
      println("total variance  " + totalVariance)

    //  val numOfSamplePoints = i
    //  println(s"num of sampled points: $numOfSamplePoints")
   //   val sampler = UniformMeshSampler3D(referenceMesh, numberOfPoints = 50000)
      val lowRankGP = LowRankGaussianProcess.approximateGPCholesky(referenceMesh.pointSet,gp,0.02, NearestNeighborInterpolator())
    //  val lowRankGP = LowRankGaussianProcess.approximateGPNystrom(gp, sampler, numBasisFunctions = 100)

      val approximatedVar = lowRankGP.klBasis.map(_.eigenvalue).sum
      println(s"Ratio of approximated variance " + approximatedVar / totalVariance)

      println(lowRankGP.klBasis.map(_.eigenvalue))
      val mm = StatisticalMeshModel(referenceMesh, lowRankGP)
      val modelGroup = ui.createGroup(s"Model-$i")
      ui.show(modelGroup, mm, "model")
      println("Number of points",mm.referenceMesh.pointSet.numberOfPoints)
      StatismoIO.writeStatismoMeshModel(mm, outputModelFile)
    }
  }
}
