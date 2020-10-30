package apps.teeth.registeration

import java.io.File

//import apps.teeth.registeration.NonRigidRegistration.fitModel
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.{LoadTestData, alignShapes}
import apps.util.AlignmentTransforms
import breeze.linalg.DenseVector
import breeze.stats.distributions.Uniform
import scalismo.common.{Domain, NearestNeighborInterpolator, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, LBFGSOptimizer, Sampler}
import scalismo.registration.{GaussianProcessTransformationSpace, L2Regularizer, MeanHuberLossMetric, MeanSquaresMetric, Registration}
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.{DiscreteLowRankGPTransformationView, ScalismoUI, ShowInScene}
import scalismo.utils.Random
import spire.syntax.rng


object GradientbasedOptimization {
  implicit val rng = scalismo.utils.Random(42)
     val ui = ScalismoUI()
  def main(args: Array[String]):Unit ={

    scalismo.initialize()

  //  val ui = ScalismoUI()

    val (model, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget("Partial")
    val alignedMeshLms = alignShapes(modelLms,targetMeshes,targetLms).align()
    val mesh = MeshIO.readMesh(new File(generalPath + s"/Registered/partialMeshes/TargetToModelFit_1.vtk")).get
    val modelGroup = ui.createGroup("ModelGroup")



    val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
    val lowRankGP = model.gp.interpolate(interpolator)
    val model2= StatisticalMeshModel.apply(mesh,lowRankGP)
    val gpView =  ui.show(modelGroup,model2, "moving mesh")

    val initialCoefficients = DenseVector.zeros[Double](model.rank)
    val registrationParameters = Seq(
      RegistrationParameters(regularizationWeight = 0.5, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.1, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.05, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.01, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.005, numberOfIterations = 1000, numberOfSampledPoints = 2000),
      RegistrationParameters(regularizationWeight = 0.001, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.0005, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.0001, numberOfIterations = 1000, numberOfSampledPoints = 1000),
      RegistrationParameters(regularizationWeight = 0.00001, numberOfIterations = 1000, numberOfSampledPoints = 1000)
    )

    alignedMeshLms.foreach{ m =>
      if( m._1._2 == "0012_36.ply") {
        val f = m._1._1

        ui.show(f, s"Fixed mesh ${m._1._2}")
        println(s"Fixed mesh ${m._1._2}")

        print("new coefficients")
        gpView.shapeModelTransformationView.shapeTransformationView.coefficients = model2.coefficients(mesh)
        val finalCoefficients: DenseVector[Double] = registrationParameters.foldLeft(initialCoefficients)((modelCoefficients, regParameters) =>
          doRegistration(model2.gp.interpolate(interpolator), mesh, f, regParameters, modelCoefficients, gpView))
      }
    }

  }

  case class RegistrationParameters(regularizationWeight : Double, numberOfIterations : Int, numberOfSampledPoints : Int)
  def doRegistration(lowRankGP : LowRankGaussianProcess[_3D, EuclideanVector[_3D]],
                     referenceMesh :TriangleMesh[_3D],
                     targetmesh : TriangleMesh[_3D],
                     registrationParameters : RegistrationParameters,
                     initialCoefficients : DenseVector[Double],
                     gpView: ShowInScene.ShowInSceneStatisticalMeshModel.View) : DenseVector[Double] =
  {
      implicit val rng = scalismo.utils.Random(42)
    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)

    val sampler1 =  PointsWithLikelyCorrespondenceSampler2(targetmesh,referenceMesh,1.0)
   //   val sampler2=  FixedPointsUniformMeshSampler3D(targetmesh,1000)
    val s = sampler1.sample()
    val pts2 =   s.map(p => referenceMesh.pointSet.findClosestPoint(p._1).id).toIndexedSeq
     val pts =   s.map(p => p._1).toIndexedSeq
   // val m = referenceMesh.marginal(pts2)
  //  ui.show(UnstructuredPointsDomain(pts),"sd").radius = 0.1
    val fixedImage =   targetmesh.operations.toDistanceImage
     val movingImage =  referenceMesh.operations.toDistanceImage

    val metric1 = MeanHuberLossMetric(
       movingImage,fixedImage,
      transformationSpace,
      sampler1
    )
    val optimizer1 = LBFGSOptimizer(registrationParameters.numberOfIterations)
    val regularizer1 = L2Regularizer(transformationSpace)
    val registration = Registration(
      metric1,
      regularizer1,
      registrationParameters.regularizationWeight,
      optimizer1
    )
    val registrationIterator = registration.iterator(initialCoefficients)
    val visualizingRegistrationIterator = for ((it, itnum) <- registrationIterator.zipWithIndex) yield {
      println(s"object value in iteration $itnum is ${it.value}")
      gpView.shapeModelTransformationView.shapeTransformationView.coefficients = it.parameters
      it
    }
    val registrationResult = visualizingRegistrationIterator.toSeq.last
    registrationResult.parameters
  }


}
case class PointsWithLikelyCorrespondenceSampler2(targetMesh: TriangleMesh[_3D],
                                                  refmesh: TriangleMesh[_3D],
                                                  maxMd: Double)     (implicit rand: Random)
  extends Sampler[_3D] {

  val ptsWithDist = refmesh.pointSet.points.toIndexedSeq.zipWithIndex
    .map {
      case (refPt,refPtId) =>
          val id =
          refmesh.pointSet.findClosestPoint(refPt).id
          val intersection = refmesh.operations.getIntersectionPoints(refPt,refmesh.vertexNormals.atPoint(id))
          val dist = intersection.minBy{ p =>
          (p - refPt).norm
          }
          (dist, 1.0)

    }

  val pts = ptsWithDist
  //  .filter {
  //    case (refPt, dist) =>  dist < maxMd
 //   }
    .map {
      case (refPt, dist) => (refPt, 1.0)
    }
    .map {
      case (refPt, dist) => (refPt, 1.0)
    }
    .toIndexedSeq
  val newpts = pts
    .filter{
      case(p,d) =>
        val bool = Math.round(targetMesh.operations.shortestDistanceToSurfaceSquared(p))
        val id =  refmesh.pointSet.findClosestPoint(p).id
        val intersection = targetMesh.operations.hasIntersection(p,refmesh.vertexNormals.atPoint(id))
        val domain = Domain.intersection(targetMesh.operations.toDistanceImage.domain, refmesh.operations.toDistanceImage.domain)

        domain.isDefinedAt(p) 

    }.map {
    case (refPt, dist) => (refPt, 1.0)
  }
  override val volumeOfSampleRegion = 1.0
  override val numberOfPoints = newpts.size
  override def sample() = {
    println(s"Sampled: $numberOfPoints");
    val pts2 =   newpts.map(p =>( refmesh.operations.closestPointOnSurface(p._1).point,1.0)).toIndexedSeq
    val distrDim1 = Uniform(0, 1000)(rand.breezeRandBasis)
    val pts = (0 until 1000).map(i => (pts2.toIndexedSeq(distrDim1.draw().toInt)))
    pts
  }
}
