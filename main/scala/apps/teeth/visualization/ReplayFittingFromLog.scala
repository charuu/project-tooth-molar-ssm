package apps.teeth.visualization

import java.awt.Color
import java.io.File

import api.sampling.ModelFittingParameters
import api.sampling.loggers.{JSONAcceptRejectLogger, jsonLogFormat}
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.LoadTestData
import scalismo.mesh.MeshMetrics
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

object ReplayFittingFromLog {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val logPath = new File(generalPath, "log")

    val (model, modelLm, targetMesh, targetLM) = LoadTestData.modelAndTarget("Partial")

    val jsonFileName = "icpProposalRegistration.json"

    println(new File(logPath, jsonFileName).toString)

    val logObj = new JSONAcceptRejectLogger[ModelFittingParameters](new File(logPath, jsonFileName))
    val logInit: IndexedSeq[jsonLogFormat] = logObj.loadLog()

    val ui = ScalismoUI(jsonFileName)
    val targetGroup = ui.createGroup("target")
    val modelGroup = ui.createGroup("model")

    ui.show(targetGroup, targetMesh.seq(0)._1, "target").color = Color.YELLOW
 ///   ui.show(modelLm,"ModelLM")
 //   ui.show(targetLM,"targetLM")
    val modelShow = ui.show(modelGroup, model, "model")

    def getLogIndex(i: Int): Int = {
      if (logInit(i).status) i
      else getLogIndex(i - 1)
    }

    val firstIndexNotReject = logInit.filter(f => f.status).head.index

    val takeEveryN = 5

    println(s"takeEvery: $takeEveryN total log : " + logInit.length)
    Thread.sleep(3000)
    for (cnt <- firstIndexNotReject until logInit.length by takeEveryN) yield {
      val index = getLogIndex(cnt)
      println(s"Index from Markov-Chain: $cnt, Index closest accepted sample: $index")
      val js = logInit(index)
      val pars = logObj.sampleToModelParameters(js)
      Thread.sleep(1000)
      val rigidTrans = ModelFittingParameters.poseTransform(pars)
      modelShow.shapeModelTransformationView.poseTransformationView.transformation = rigidTrans
      modelShow.shapeModelTransformationView.shapeTransformationView.coefficients = pars.shapeParameters.parameters
      Thread.sleep(1000)
      println(MeshMetrics.avgDistance(targetMesh.seq(0)._1, model.instance(pars.shapeParameters.parameters).transform(rigidTrans)))
    }
  }
}
