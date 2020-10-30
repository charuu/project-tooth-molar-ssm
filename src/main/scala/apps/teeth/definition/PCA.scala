package apps.teeth.definition

import java.io.File

import apps.teeth.utilities.Paths.generalPath
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.ui.api.ScalismoUI

object PCA {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    implicit val rng = scalismo.utils.Random(42)
    val ui = ScalismoUI()
    val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Registered/model/AugmentedPCAD.h5")).get
   // val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Ref/model/tooth_gp_model_A-components.h5")).get

    val registeredMeshesFileList = new java.io.File(generalPath + "/Registered/partialMeshes/ICPProposal/").listFiles
   // val registeredMeshesFileList = new java.io.File(generalPath + "/Registered/partialMeshes/PCA/").listFiles

    val registeredMeshes = (registeredMeshesFileList).map(m => MeshIO.readMesh(m).get)

    registeredMeshesFileList.map{ mesh =>
      println(MeshIO.readMesh(mesh).get.pointSet.numberOfPoints,": number of points ")
      println(model.referenceMesh.pointSet.numberOfPoints,": number of points ")
      ui.show(MeshIO.readMesh(mesh).get,s"mesh-${mesh.getName}")
    }
    val modelGroup1 = ui.createGroup("modelGroup")
    val dc = DataCollection.fromMeshSequence(model.referenceMesh, registeredMeshes)._1.get
    val alignedDC = DataCollection.gpa(dc)

    val pcaModel = StatisticalMeshModel.createUsingPCA(alignedDC).get
    ui.show(modelGroup1, pcaModel, "PcaModel")

    StatisticalModelIO.writeStatisticalMeshModel(pcaModel,new File(generalPath + "/Registered/model/ssModel.h5"))

  }
  }
