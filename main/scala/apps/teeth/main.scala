package apps.teeth

import apps.teeth.utilities.{LoadTestData, alignShapes}
import scalismo.ui.api.ScalismoUI

object main {
  def main(args: Array[String]):Unit = {
    scalismo.initialize()
    implicit val rng = scalismo.utils.Random(42)
    val ui = ScalismoUI(s"ScalismoUI")
    //load full samples data
    val (model, modelLms, sampleMeshes, sampleLMs) = LoadTestData.modelAndTarget("Full")
    //align meshes
    val Seq(alignedMesh,alignedLm) = alignShapes(modelLms,sampleMeshes,sampleLMs).align()

    //registeration

    //--non rigid icp
    //--grad opt
    //--mcmc
    //build pca model
    //Augment pca model
    //load partial data
    //partial registeration
    //--non rigid icp
    //-- icp proposal mcmc
    //pca model

  }
}
