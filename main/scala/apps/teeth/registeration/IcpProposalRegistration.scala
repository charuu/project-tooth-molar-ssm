package apps.teeth.registeration

import java.awt.Color
import java.io.File

import api.other.{RegistrationComparison, TargetSampling}
import api.sampling.evaluators.TargetToModelEvaluation
import api.sampling.{MixedProposalDistributions, ModelFittingParameters, ProductEvaluators, SamplingRegistration}
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.LoadTestData
import scalismo.geometry._3D
import scalismo.mesh.{MeshMetrics, TriangleMesh, TriangleMesh3D}
import scalismo.sampling.DistributionEvaluator
import scalismo.sampling.proposals.MixtureProposal.ProposalGeneratorWithTransition
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless, StatisticalMeshModelViewControls}

object IcpProposalRegistration {

  def fitting(model: StatisticalMeshModel, targetMesh: TriangleMesh3D, evaluator: Map[String, DistributionEvaluator[ModelFittingParameters]], proposal: ProposalGeneratorWithTransition[ModelFittingParameters], numOfIterations: Int, showModel: Option[StatisticalMeshModelViewControls], log: File, initialParameters: Option[ModelFittingParameters] = None): TriangleMesh[_3D] = {

    val samplingRegistration = new SamplingRegistration(model, targetMesh, showModel, modelUiUpdateInterval = 5, acceptInfoPrintInterval = 5)
    val t0 = System.currentTimeMillis()

    val best = samplingRegistration.runfitting(evaluator, proposal, numOfIterations, initialModelParameters = initialParameters, jsonName = log)

    val t1 = System.currentTimeMillis()
    println(s"ICP-Timing: ${(t1 - t0) / 1000.0} sec")
    ModelFittingParameters.transformedMesh(model, best)
  }

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUIHeadless()

    println(s"Starting Metropolis Hastings registrations with ICP-proposal!")

    val logPath = new File(generalPath, "log" +
      "")

    val (m, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget("Partial")
   // model.truncate(50)
    val targetMesh = targetMeshes.seq(0)._1.operations.decimate(5000)
    val model = m.decimate(5000)
    val numOfEvaluatorPoints = 500 // Used for the likelihood evaluator
    val numOfICPPointSamples = 500 // Used for the ICP proposal
    val numOfSamples = 5000 // Length of Markov Chain

  //  val proposalIcp = MixedProposalDistributions.mixedProposalRandom(model)
   val proposalIcp = MixedProposalDistributions.mixedProposalICPGenerator(model, targetMesh, numOfICPPointSamples, projectionDirection = TargetSampling,0.2,0.2,0.1)

    // Euclidean likelihood evaluator using a Gaussian distribution
    val evaluator = ProductEvaluators.proximityAndIndependent(model, targetMesh, TargetToModelEvaluation, uncertainty = 2.0, numberOfEvaluationPoints = numOfEvaluatorPoints)
    // Hausdorff likelihood evaluator using an Exponential distribution
    //    val evaluator = ProductEvaluators.proximityAndHausdorff(model, targetMesh, uncertainty = 100.0)


    val modelGroup = ui.createGroup("modelGroup")
    val targetGroup = ui.createGroup("targetGroup")
    val finalGroup = ui.createGroup("finalGroup")

    val showModel = ui.show(modelGroup, model, "model")
    ui.show(modelGroup, modelLms, "landmarks")
    val showTarget = ui.show(targetGroup, targetMesh, "target")
    ui.show(targetGroup, targetLms.seq(0)._1, "landmarks")
    showTarget.color = Color.YELLOW


    val bestRegistration = fitting(model, targetMesh, evaluator, proposalIcp, numOfSamples, Option(showModel), new File(logPath, s"icpProposalRegistration.json"))
    ui.show(finalGroup, bestRegistration, "best-fit")
    println(MeshMetrics.avgDistance(targetMesh,bestRegistration))
    RegistrationComparison.evaluateReconstruction2GroundTruth("SAMPLE", bestRegistration ,targetMesh)
  }
}
