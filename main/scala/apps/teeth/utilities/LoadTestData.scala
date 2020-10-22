package apps.teeth.utilities

import java.io.File

import Paths.generalPath
import apps.util.AlignmentTransforms
import scalismo.geometry.{EuclideanVector, Landmark, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatismoIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.utils.Random

object LoadTestData {
  implicit val random: Random = Random(1024)

  def modelAndTarget(meshType: String): (StatisticalMeshModel, Seq[(Landmark[_3D])], Seq[(TriangleMesh[_3D],String)], Seq[(Seq[Landmark[_3D]],String)]) = {

    val (modelFile, modelLmsFile, targetMeshFile, targetLmsFile) = LoadFiles(meshType)
    val m = StatismoIO.readStatismoMeshModel(modelFile).get
    val model = m.truncate(300)
    val modelLms = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get

    println(s"Model file to be used: $modelFile")
    val sampleMeshes = targetMeshFile.map(m => (MeshIO.readMesh(m).get, m.getName()))
    val sampleMeshesLandmarks = targetLmsFile.map(l => (LandmarkIO.readLandmarksJson[_3D](l).get, l.getName()))
    (model, modelLms, sampleMeshes, sampleMeshesLandmarks)

  }

  def LoadFiles(meshType: String): (File,File,Seq[File],Seq[File]) = meshType match {
    case "Partial"=>
      val modelFile = new File(generalPath, "Registered/model/AugmentedPCA_0606.h5")
      val modelLmsFile = new File(generalPath, "Registered/model/AugmentedPCA.json")
      val targetMeshFile = new File(generalPath, "Raw/partial/mirrored_meshes/").listFiles().toSeq
      val targetLmsFile = new File(generalPath, s"Raw/partial/landmarks/").listFiles().toSeq
      (modelFile, modelLmsFile, targetMeshFile, targetLmsFile)
    case "Full" =>
      val modelFile = new File(generalPath, "/Ref/model/tooth_gp_model_50000-components.h5")
      val modelLmsFile = new File(generalPath, "/Ref/landmarks/ref_landmarks_7.json")
      val sampleMeshFileList = new File(generalPath, "/Raw/mesh/10_Sample_meshes/Meshes/").listFiles().toSeq
      val sampleMeshLandmarksFileList = new File(generalPath, "/Raw/landmarks/noisyMeshLandmarks/").listFiles().toSeq
      (modelFile, modelLmsFile, sampleMeshFileList, sampleMeshLandmarksFileList)

  }
}

