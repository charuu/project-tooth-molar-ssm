
import java.io.File

import apps.teeth.utilities.Paths.generalPath
import apps.teeth.registeration.RegistrationFullMesh.model
import apps.teeth.utilities.LoadTestData
import breeze.linalg.DenseVector
import scalismo.common.{Field, NearestNeighborInterpolator, RealSpace}
import scalismo.geometry.{EuclideanVector, IntVector, NDSpace, Point, _3D}
import scalismo.image.{DifferentiableScalarImage, DiscreteImageDomain, ScalarImage}
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.mesh.TriangleMesh
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, GridSampler, LBFGSOptimizer, PointsWithLikelyCorrespondenceSampler, RandomMeshSampler3D, Sampler, UniformMeshSampler3D}
import scalismo.registration.{GaussianProcessTransformationSpace, L2Regularizer, MeanHuberLossMetric, MeanPointwiseLossMetric, MeanSquaresMetric, Registration, TransformationSpace}
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI

object GradientbasedOptimizationForPartial {
  implicit val rng = scalismo.utils.Random(42)
  scalismo.initialize()

  val ui = ScalismoUI()
  val (model, modelLms, targetMesh, targetLms) = LoadTestData.modelAndTarget2()
  val modelFile = new File(generalPath, "Registered/model/AugmentedPCA011000102001000.h5")

  val modelGroup = ui.createGroup("model")

  val gpView = ui.addTransformation(modelGroup, model.gp, "gp")
  val covSSM : MatrixValuedPDKernel[_3D] =  DiagonalKernel[_3D](GaussianKernel(2), 3) * 1.0
  val zeroMean = Field(RealSpace[_3D], (pt: Point[_3D]) => EuclideanVector(0, 0, 0))

  val gp = GaussianProcess(zeroMean, covSSM)

  def main(args: Array[String]):Unit ={
    val alignedMeshes = targetMesh
    val sampleMeshesFileList = new java.io.File(generalPath+"/aligned/meshes").listFiles.sortBy(f => f.getName())

    ui.show(alignedMeshes, "Fixed mesh")

    val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
    val lowRankGP = model.gp.interpolate(interpolator)


    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)

    ui.show(modelGroup,model, "moving mesh")


    val initialCoefficients = DenseVector.zeros[Double](model.rank)
    val registrationParameters = Seq(
      RegistrationParameters(regularizationWeight = 4.0, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 2.0, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 1.0, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.01, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.001, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.0001, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.00001, numberOfIterations = 500, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.000001, numberOfIterations = 500, numberOfSampledPoints = 1000)
    )
   // sampleMeshesFileList.map{ m =>
      val finalCoefficients = registrationParameters.foldLeft(initialCoefficients)((modelCoefficients, regParameters) =>
        doRegistration(lowRankGP, model.referenceMesh,alignedMeshes, regParameters, modelCoefficients))

      val registrationTransformation = transformationSpace.transformForParameters(finalCoefficients)
      val targetMeshOperations = alignedMeshes
      val projection = (pt : Point[_3D]) => {
        targetMeshOperations.pointSet.findClosestPoint(pt).point
      }

      val finalTransformation = registrationTransformation.andThen(projection)

      val projectedMesh = model.referenceMesh.transform(finalTransformation)
      val resultGroup = ui.createGroup("result")
      val projectionView = ui.show(resultGroup, projectedMesh, "projection")
      MeshIO.writeSTL(projectedMesh, new File(generalPath +s"/Registered/mesh/0012_36.stl"))

  //  }

  }

  case class RegistrationParameters(regularizationWeight : Double, numberOfIterations : Int, numberOfSampledPoints : Int)
  def doRegistration(
                      lowRankGP : LowRankGaussianProcess[_3D, EuclideanVector[_3D]],
                      referenceMesh : TriangleMesh[_3D],
                      targetmesh : TriangleMesh[_3D],
                      registrationParameters : RegistrationParameters,
                      initialCoefficients : DenseVector[Double]
                    ) : DenseVector[Double] =
  {
    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)
    val fixedImage = referenceMesh.operations.toDistanceImage
    val movingImage = targetmesh.operations.toDistanceImage

    val sampler1 = PointsWithLikelyCorrespondenceSampler2(
      gp,referenceMesh,targetmesh,2.0
    )

    val metric = MeanHuberLossMetric(
      fixedImage,
      movingImage,
      transformationSpace,
      sampler1
    )
    val optimizer = LBFGSOptimizer(registrationParameters.numberOfIterations)
    val regularizer = L2Regularizer(transformationSpace)
    val registration = Registration(
      metric,
      regularizer,
      registrationParameters.regularizationWeight,
      optimizer
    )
    val registrationIterator = registration.iterator(initialCoefficients)
    val visualizingRegistrationIterator = for ((it, itnum) <- registrationIterator.zipWithIndex) yield {
      println(s"object value in iteration $itnum is ${it.value}")
      gpView.coefficients = it.parameters
      it
    }
    val registrationResult = visualizingRegistrationIterator.toSeq.last
    registrationResult.parameters
  }

}
case class PointsWithLikelyCorrespondenceSampler2(gp: GaussianProcess[_3D, EuclideanVector[_3D]], refmesh: TriangleMesh[_3D], targetMesh: TriangleMesh[_3D], maxMd: Double) extends Sampler[_3D] {

  //  val meanPts = refmesh.points.map(gp.mean(_).toPoint)
  val meanPts = refmesh.pointSet.points.map {
    x: Point[_3D] => x + gp.mean(x)
  }
  val ptsWithDist = refmesh.pointSet.points.toIndexedSeq.zipWithIndex.par
    .map {
      case (refPt, refPtId) =>

        val closestTgtPt = targetMesh.pointSet.findClosestPoint(refPt).point
        (refPt, gp.marginal(refPt).mahalanobisDistance((closestTgtPt - refPt).toBreezeVector))
    }

  val pts = ptsWithDist
    .filter {
      case (refPt, dist) => dist < maxMd
    }
    .map {
      case (refPt, dist) => (refPt, 1.0)
    }
    .map {
      case (refPt, dist) => (refPt, 1.0)
    }
    .toIndexedSeq

  override val volumeOfSampleRegion = 1.0
  override val numberOfPoints = pts.size
  override def sample() = {
    println(s"Sampled: $numberOfPoints"); pts
  }
}

