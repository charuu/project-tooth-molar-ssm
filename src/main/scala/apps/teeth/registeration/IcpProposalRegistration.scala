package apps.teeth.registeration

import java.awt.Color
import java.io.File

import api.other.{DoubleProjectionSampling, ModelAndTargetSampling, ModelSampling, RegistrationComparison, TargetSampling}
import api.sampling.evaluators.{SymmetricEvaluation, TargetToModelEvaluation}
import api.sampling.{MixedProposalDistributions, ModelFittingParameters, ProductEvaluators, SamplingRegistration, SurfaceNoiseHelpers}
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.{LoadTestData, alignShapes}
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
    val meshType = "Full"
    println(s"Starting Metropolis Hastings registrations with ICP-proposal!")

    val (m, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget(meshType)
    val alignedMeshLms = alignShapes(modelLms, targetMeshes, targetLms).align()
    val logPath = new File(generalPath, "log" + "")

    alignedMeshLms.map { mesh =>
        println("Targetmesh :", mesh._1._2)
        val targetMesh = mesh._1._1.operations.decimate(10000)

        val model = m.decimate(10000)
        val numOfEvaluatorPoints = model.referenceMesh.pointSet.numberOfPoints // Used for the likelihood evaluator
        val numOfICPPointSamples = numOfEvaluatorPoints // Used for the ICP proposal
        val numOfSamples = 300 // Length of Markov Chain

        val proposalIcp =     MixedProposalDistributions.mixedProposalICPGenerator(model, targetMesh, numOfICPPointSamples, projectionDirection = ModelAndTargetSampling, 0.2, 0.2, 0.1)

        val evaluator = ProductEvaluators.proximityAndIndependent(model, targetMesh, SymmetricEvaluation, uncertainty = 0.0001, numberOfEvaluationPoints = numOfEvaluatorPoints, "Surface")

      //utils.Plot.plotTrace(samples,evaluator,new java.io.File("data/handedData/plotTraceTotal.png"))
        val modelGroup = ui.createGroup("modelGroup")
        val targetGroup = ui.createGroup("targetGroup")
        val finalGroup = ui.createGroup("finalGroup")

        val showModel = ui.show(modelGroup, model, "model")
       // ui.show(modelGroup, modelLms, "landmarks")
        val showTarget = ui.show(targetGroup, targetMesh, "target")
      //  ui.show(targetGroup, alignedMeshLms.seq(0)._2._1, "landmarks")
        showTarget.color = Color.YELLOW

        val bestRegistration = fitting(model, targetMesh, evaluator, proposalIcp, numOfSamples, Option(showModel), new File(logPath, s"icpProposalRegistration ${mesh._1._2.replace(".stl", "")}.json"))

     //   ui.show(finalGroup, bestRegistration, "best-fit")
        println(MeshMetrics.avgDistance(targetMesh, bestRegistration))
        RegistrationComparison.evaluateReconstruction2GroundTruth("SAMPLE", bestRegistration, targetMesh)
      }
    }
 // }
}
