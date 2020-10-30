package apps.teeth.visualization

import java.awt.Color
import java.io.File

import api.sampling.{ModelFittingParameters, SurfaceNoiseHelpers}
import api.sampling.loggers.{JSONAcceptRejectLogger, jsonLogFormat}
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.{LoadTestData, alignShapes}
import scalismo.mesh.MeshMetrics
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

object ReplayFittingFromLog {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val logPath = new File(generalPath, "log")

    val (m, modelLm, targetMesh, targetLM) = LoadTestData.modelAndTarget("Partial")
    val alignedMeshLms = alignShapes(modelLm, targetMesh, targetLM).align()

    alignedMeshLms.map { target=>
      if (target._1._2 == "0056_36.ply") {

        println("Targetmesh :", target._1._2)

      //  val jsonFileName = s"icpProposalRegistration ${target._1._2.replace(".stl",".json")}"
        val jsonFileName = s"icpProposalRegistration 0056_36.json"

        println(new File(logPath, jsonFileName).toString)

        val logObj = new JSONAcceptRejectLogger[ModelFittingParameters](new File(logPath, jsonFileName))
        val logInit: IndexedSeq[jsonLogFormat] = logObj.loadLog()

        val ui = ScalismoUI(jsonFileName)
        val targetGroup = ui.createGroup("target")
        val modelGroup = ui.createGroup("model")

        ui.show(targetGroup, target._1._1, "target").color = Color.YELLOW

        val modelShow = ui.show(modelGroup, m, "model")

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
          println(MeshMetrics.avgDistance(m.instance(pars.shapeParameters.parameters).transform(rigidTrans), target._1._1))
        }
      }
    }
  }
}
