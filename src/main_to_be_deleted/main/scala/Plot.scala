/*
package utils

import java.io.File

import api.sampling.ModelFittingParameters
import com.cibo.evilplot.numeric.Point
import com.cibo.evilplot.plot.{BoxPlot, LinePlot}
import com.cibo.evilplot.plot.aesthetics.DefaultTheme._
import scalismo.sampling.evaluators.ProductEvaluator

import scala.collection.immutable

object Plot {
  def plotParameters(samplesASM: IndexedSeq[ModelFittingParameters], file: File) = {
    val numCoeffs = samplesASM.head.shapeParameters.parameters.data.length
    val coeffs: Seq[Seq[Double]] = for (i <- 0 until numCoeffs) yield {
      for (s <- samplesASM) yield {
        s.shapeParameters.parameters.data(i)
      }
    }


    val plot = BoxPlot(coeffs)
      .standard(xLabels = (1 to numCoeffs).map(_.toString))
      .ybounds(-3.0, 3.0)
      .render()
      .asBufferedImage
    javax.imageio.ImageIO.write(plot, "png", file)

  }

  def plotTrace(samplesASM: IndexedSeq[ModelFittingParameters], asmPosteriorEvaluator: ProductEvaluator[ModelFittingParameters], file: java.io.File) = {

    val logValue = samplesASM.map(asmPosteriorEvaluator.logValue)

    val plot = LinePlot(logValue.zipWithIndex.map { case (v, i) => Point(i.toDouble, v) })
      .render()
      .asBufferedImage
    javax.imageio.ImageIO.write(plot, "png", file)
  }
}*/
