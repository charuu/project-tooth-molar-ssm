package apps.teeth.utilities

import api.sampling.ModelFittingParameters
import scalismo.sampling.loggers.AcceptRejectLogger
import scalismo.sampling.{DistributionEvaluator, ProposalGenerator}

class Logger extends AcceptRejectLogger[ModelFittingParameters] {
  private val numAccepted = collection.mutable.Map[String, Int]()
  private val numRejected = collection.mutable.Map[String, Int]()

  override def accept(current: ModelFittingParameters,
                      sample: ModelFittingParameters,
                      generator: ProposalGenerator[ModelFittingParameters],
                      evaluator: DistributionEvaluator[ModelFittingParameters]
                     ): Unit = {
    val numAcceptedSoFar = numAccepted.getOrElseUpdate(sample.generatedBy, 0)
    println("Accepted Sample: ",evaluator.logValue(sample)," ,current : ",evaluator.logValue(current))
    numAccepted.update(sample.generatedBy, numAcceptedSoFar + 1)
  }

  override def reject(current: ModelFittingParameters,
                      sample: ModelFittingParameters,
                      generator: ProposalGenerator[ModelFittingParameters],
                      evaluator: DistributionEvaluator[ModelFittingParameters]
                     ): Unit = {
    val numRejectedSoFar = numRejected.getOrElseUpdate(sample.generatedBy, 0)
    println("Rejected Sample: ",evaluator.logValue(sample)," ,current : ",evaluator.logValue(current))
    numRejected.update(sample.generatedBy, numRejectedSoFar + 1)
  }


  def acceptanceRatios() : Map[String, Double] = {
    val generatorNames = numRejected.keys.toSet.union(numAccepted.keys.toSet)
    val acceptanceRatios = for (generatorName <- generatorNames ) yield {
      val total = (numAccepted.getOrElse(generatorName, 0)
        + numRejected.getOrElse(generatorName, 0)).toDouble
      (generatorName, numAccepted.getOrElse(generatorName, 0) / total)
    }
    acceptanceRatios.toMap
  }
}