package apps.teeth.visualization

import java.io.File

import api.sampling.ModelFittingParameters
import api.sampling.loggers.{JSONAcceptRejectLogger, jsonLogFormat}
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.LoadTestData
import apps.util.{LogHelper, PosteriorVariability}
import scalismo.ui.api.ScalismoUI

object PosteriorVariabilityToMeshColor {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val logPath = new File(generalPath, "log")

    val (model, _, _, _) = LoadTestData.modelAndTarget("Partial")

    val jsonFileName = "icpProposalRegistration.json"

    println(new File(logPath, jsonFileName).toString)

    val logObj = new JSONAcceptRejectLogger[ModelFittingParameters](new File(logPath, jsonFileName))
    val logInit: IndexedSeq[jsonLogFormat] = logObj.loadLog()
    val burnInPhase = 10

    val logSamples = LogHelper.samplesFromLog(logInit, takeEveryN = 50, total = 10000, burnInPhase)
    println(s"Number of samples from log: ${logSamples.length}/${logInit.length - burnInPhase}")
    val logShapes = LogHelper.logSamples2shapes(model, logSamples.map(_._1))

    val best = ModelFittingParameters.transformedMesh(model, logObj.getBestFittingParsFromJSON)

    val colorMap_normalVariance = PosteriorVariability.computeDistanceMapFromMeshesNormal(logShapes, best, sumNormals = true)
    val colorMap_posteriorEstimate = PosteriorVariability.computeDistanceMapFromMeshesTotal(logShapes, best)

    val ui = ScalismoUI(s"Posterior visualization - $jsonFileName")
    val modelGroup = ui.createGroup("model")
    val colorGroup = ui.createGroup("color")
    val showModel = ui.show(modelGroup, model, "model")
    showModel.meshView.opacity = 0.0
    ui.show(colorGroup, colorMap_posteriorEstimate, "posterior")
    ui.show(colorGroup, colorMap_normalVariance, "normal")
    ui.show(colorGroup, best, "best-fit")


  }
}
