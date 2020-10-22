import java.io.File

import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.ui.api.ScalismoUI
import main.scala.Paths.generalPath
import scalismo.common.NearestNeighborInterpolator
import scalismo.geometry._3D
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel}
import scalismo.numerics.UniformMeshSampler3D
object PCA {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    implicit val rng = scalismo.utils.Random(42)
    val ui = ScalismoUI()
    val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Ref/model/tooth_gp_model_1000-components.h5")).get

    val registeredMeshesFileList = new java.io.File(generalPath + "/Registered/mesh/1000_points_registered").listFiles
    val registeredMeshes = registeredMeshesFileList.map(m => MeshIO.readMesh(m).get)

    registeredMeshesFileList.map{ mesh =>
      println(MeshIO.readMesh(mesh).get.pointSet. ,": number of points ")
      ui.show(MeshIO.readMesh(mesh).get,s"mesh-${mesh.getName}")
    }
    val modelGroup1 = ui.createGroup("modelGroup")
    val dc = DataCollection.fromMeshSequence(model.referenceMesh, registeredMeshes)._1.get
    val pcaModel = StatisticalMeshModel.createUsingPCA(dc).get
    ui.show(modelGroup1, pcaModel, "PcaModel")

    StatisticalModelIO.writeStatisticalMeshModel(pcaModel,new File(generalPath + "/Registered/model/pcaModel_1000points_1205.h5"))

  }
  }