
import java.io.File

import apps.teeth.utilities.Paths.generalPath
import apps.teeth.registeration.RegistrationFullMesh.model
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
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}

object GradientbasedOptimization {
 implicit val rng = scalismo.utils.Random(42)
 scalismo.initialize()

 val ui = ScalismoUI()

 val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath +"/Ref/model/tooth_gp_model_1000-components.h5")).get
 val modelGroup = ui.createGroup("model")

 val gpView = ui.addTransformation(modelGroup, model.gp, "gp")

 def main(args: Array[String]):Unit ={
  val sampleMeshesFileList = new java.io.File(generalPath+"/aligned/meshes").listFiles.sortBy(f => f.getName())


  val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
  val lowRankGP = model.gp.interpolate(interpolator)
  val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)

  ui.show(modelGroup,model, "moving mesh")


  val initialCoefficients = DenseVector.zeros[Double](model.rank)
  val registrationParameters = Seq(
   RegistrationParameters(regularizationWeight = 2.0, numberOfIterations = 20, numberOfSampledPoints = 1000),
   RegistrationParameters(regularizationWeight = 1e-1, numberOfIterations = 20, numberOfSampledPoints = 1000),
   RegistrationParameters(regularizationWeight = 1e-2, numberOfIterations = 100, numberOfSampledPoints = 1000),
   RegistrationParameters(regularizationWeight = 1e-4, numberOfIterations = 200, numberOfSampledPoints = 1000),
   RegistrationParameters(regularizationWeight = 1e-6, numberOfIterations = 1000, numberOfSampledPoints = 1000),
   RegistrationParameters(regularizationWeight = 1e-8, numberOfIterations = 1000, numberOfSampledPoints = 1000)
  )
  sampleMeshesFileList.take(1).map{ m =>
    val fixedMesh = MeshIO.readMesh(m).get
   ui.show(fixedMesh, s"Fixed mesh ${m.getName}")

   val finalCoefficients = registrationParameters.foldLeft(initialCoefficients)((modelCoefficients, regParameters) =>
    doRegistration(lowRankGP, model.referenceMesh,fixedMesh, regParameters, modelCoefficients))

   val registrationTransformation = transformationSpace.transformForParameters(finalCoefficients)
   val targetMeshOperations = fixedMesh
   val projection = (pt : Point[_3D]) => {
    targetMeshOperations.pointSet.findClosestPoint(pt).point
   }

   val finalTransformation = registrationTransformation.andThen(projection)

   val projectedMesh = model.referenceMesh.transform(finalTransformation)
   val resultGroup = ui.createGroup("result")
   val projectionView = ui.show(resultGroup, projectedMesh, "projection")
   MeshIO.writeSTL(projectedMesh, new File(generalPath +s"/Registered/mesh/${m.getName}"))

  }

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
  val sampler1 = FixedPointsUniformMeshSampler3D(
   referenceMesh,
   registrationParameters.numberOfSampledPoints
  )


  val metric = MeanSquaresMetric(
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

