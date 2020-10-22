package main.scala

import java.io.File

import apps.teeth.utilities.LoadTestData
import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh}
import scalismo.sampling.algorithms.MetropolisHastings
import scalismo.sampling.evaluators.ProductEvaluator
import scalismo.sampling.proposals.MixtureProposal
import scalismo.sampling.{DistributionEvaluator, ProposalGenerator, TransitionProbability}
import scalismo.statisticalmodel.{MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Memoize

import scala.util.Random

object MCMC {


  def computeCenterOfMass(mesh : TriangleMesh[_3D]) : Point[_3D] = {
    val normFactor = 1.0 / mesh.pointSet.numberOfPoints
    mesh.pointSet.points.foldLeft(Point(0, 0, 0))((sum, point) => sum + point.toVector * normFactor)
  }


  case class CachedEvaluator[A](evaluator: DistributionEvaluator[A]) extends DistributionEvaluator[A] {
    val memoizedLogValue = Memoize(evaluator.logValue, 10)

     def logValue(sample: A): Double = {
      memoizedLogValue(sample)
    }
  }


  case class PriorEvaluator(model: StatisticalMeshModel)
    extends DistributionEvaluator[ModelFittingParameters] {

    val translationPrior = breeze.stats.distributions.Gaussian(0.0, 1.0)
    val rotationPrior = breeze.stats.distributions.Gaussian(0, 1.0)

    override def logValue(sample: ModelFittingParameters): Double = {
      model.gp.logpdf(sample.shapeParameters.parameters) +
        translationPrior.logPdf(sample.poseParameters.translation.x) +
        translationPrior.logPdf(sample.poseParameters.translation.y) +
        translationPrior.logPdf(sample.poseParameters.translation.z) +
        rotationPrior.logPdf(sample.poseParameters.rotation._1) +
        rotationPrior.logPdf(sample.poseParameters.rotation._2) +
        rotationPrior.logPdf(sample.poseParameters.rotation._3)
    }
  }
  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val ui = ScalismoUI()
    implicit val rng = scalismo.utils.Random(42)
    val (model, modelLms, targetMesh, targetLms) = LoadTestData.modelAndTarget2()

    val modelGroup = ui.createGroup("model")
    val modelView = ui.show(modelGroup, model, "model")
    val modelLmView = ui.show(modelGroup, modelLms, "model")
    modelView.meshView.opacity = 0.5

    val targetGroup = ui.createGroup("target")
    ui.show(targetGroup,targetMesh,"target Mesh")
    ui.show(targetGroup,targetLms,"target Mesh")

    val likelihoodIndependent = breeze.stats.distributions.Gaussian(0,2.0)

    val likelihoodEvaluator = CachedEvaluator(IndependentPointDistanceEvaluator(model,targetMesh,likelihoodIndependent,"TargetToModelEvaluation", 100))
    val priorEvaluator = CachedEvaluator(PriorEvaluator(model))

    val posteriorEvaluator = ProductEvaluator(priorEvaluator, likelihoodEvaluator)

    val nonRigidProposal = NonRigidIcpProposal(model,targetMesh,0.0,100,"TargetSampling")
    val generator = MixtureProposal.fromProposalsWithTransition(
    (1,nonRigidProposal)
    )

    val scaleParameter = ScaleParameter(0.1)
    val poseParameter = PoseParameters(EuclideanVector(0, 0, 0),(0.0, 0.0, 0.0),computeCenterOfMass(model.mean))
    val shapeParameter = ShapeParameters( DenseVector.zeros[Double](model.rank))

    val initialSample = ModelFittingParameters(scaleParameter, poseParameter,shapeParameter,"initial")
    val chain = MetropolisHastings(generator, posteriorEvaluator)
    val logger = new Logger()
    val mhIterator = chain.iterator(initialSample,logger)

    val samplingIterator = for((sample, iteration) <- mhIterator.zipWithIndex) yield {

      if (iteration % 5 == 0) {
        println("iteration " + iteration)
        modelView.shapeModelTransformationView.shapeTransformationView.coefficients = sample.shapeParameters.parameters
    //    modelView.shapeModelTransformationView.poseTransformationView.transformation = ModelFittingParameters.poseTransform(sample)
      }

      sample
    }

    val samples = samplingIterator.take(2000).toIndexedSeq
   // utils.Plot.plotParameters(samples,new java.io.File("data/handedData/plotparametersTotal.png"))
   // utils.Plot.plotTrace(samples,posteriorEvaluator,new java.io.File("data/handedData/plotTraceTotal.png"))
    graph(samples, priorEvaluator, priorEvaluator, likelihoodEvaluator, chain.evaluator, 2)

    println("acceptance ratio is " +logger.acceptanceRatios())
    val bestSample = samples.maxBy(posteriorEvaluator.logValue)
    val bestFit = model.instance(bestSample.shapeParameters.parameters)
    println(MeshMetrics.avgDistance(targetMesh,bestFit))
    val resultGroup = ui.createGroup("result")
    ui.show(resultGroup, bestFit, "best fit")

  }
  def graph(samples: Seq[ModelFittingParameters], posteriorEvaluator: DistributionEvaluator[ModelFittingParameters],
            priorEvaluator: DistributionEvaluator[ModelFittingParameters], likelihoodEvaluator: DistributionEvaluator[ModelFittingParameters], chaindistribution: DistributionEvaluator[ModelFittingParameters], i: Int): Unit = {
    val posterior = samples.map { s =>
      posteriorEvaluator.logValue(s)
    }


    val distribution = samples.map { s =>
      chaindistribution.logValue(s)
    }

    val prior = samples.map { s =>
      priorEvaluator.logValue(s)
    }
    val liklihood = samples.map { s =>
      likelihoodEvaluator.logValue(s)
    }


    /*val value = Histogram(distribution)
      .xAxis()
      .yAxis()
      .xLabel("Distribution PDF")
      .yLabel("Samples")
      .frame()
      .render()
    val value2 = Histogram(posterior)
      .xAxis()
      .yAxis()
      .xLabel("Posterior PDF")
      .yLabel("Samples").frame()
      .render()
    val value3 = Histogram(prior)
      .xAxis()
      .yAxis()
      .xLabel("Prior PDF")
      .yLabel("Samples")
      .frame()
      .render()
    val value4 = Histogram(liklihood)
      .xAxis()
      .yAxis()
      .xLabel("Likelihood PDF")
      .yLabel("Samples").frame()
      .render()


    javax.imageio.ImageIO.write(value.asBufferedImage, "png", new java.io.File("data/handedData/distribution.png"))
    javax.imageio.ImageIO.write(value2.asBufferedImage, "png", new java.io.File("data/handedData/posterior.png"))
    javax.imageio.ImageIO.write(value3.asBufferedImage, "png", new java.io.File("data/handedData/prior.png"))
    javax.imageio.ImageIO.write(value4.asBufferedImage, "png", new java.io.File("data/handedData/liklihood.png"))
*/
  }
}

