package apps.teeth.utilities

import java.io.File

import Paths.generalPath
import apps.util.AlignmentTransforms
import scalismo.geometry.{EuclideanVector, Landmark, Point, Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatismoIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

object LoadTestData {
  implicit val random: Random = Random(1024)
  val meshType = "Partial"
  def modelAndTarget(meshType: String): (StatisticalMeshModel, Seq[(Landmark[_3D])], Seq[(TriangleMesh[_3D],String)], Seq[(Seq[Landmark[_3D]],String)]) = {

    val (modelFile, modelLmsFile, targetMeshFile, targetLmsFile) = LoadFiles(meshType)
    val m = StatismoIO.readStatismoMeshModel(modelFile).get

    val modelLms = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get

    println(s"Model file to be used: $modelFile")
    val sampleMeshes = targetMeshFile.map(m => (MeshIO.readMesh(m).get, m.getName()))
    val sampleMeshesLandmarks = targetLmsFile.map(l => (LandmarkIO.readLandmarksJson[_3D](l).get, l.getName()))
    (m, modelLms, sampleMeshes, sampleMeshesLandmarks)

  }

  def LoadFiles(meshType: String): (File,File,Seq[File],Seq[File]) = meshType match {
    case "Partial"=>
     val modelFile = new File(generalPath, "Registered/model/AugmentedPCAD.h5")
      val modelLmsFile = new File(generalPath, "Registered/model/AugmentedPCA.json")

      val targetMeshFile = new File(generalPath, "Raw/partialMeshes/mirroredMeshes/").listFiles().toSeq
      val targetLmsFile = new File(generalPath, "Raw/partialMeshes/landmarks/").listFiles().toSeq
      (modelFile, modelLmsFile, targetMeshFile, targetLmsFile)
    case "Full" =>
      val modelFile = new File(generalPath, "/Ref/model/tooth_gp_model_A-components.h5")
      val modelLmsFile = new File(generalPath, "/Ref/landmarks/ref_landmarks.json")
    //  val sampleMeshFileList = new File(generalPath, "/Registered/fullMeshes/ICPProposal/").listFiles().toSeq
      val sampleMeshFileList = new File(generalPath, "/Raw/fullMeshes/meshes/").listFiles().toSeq
      val sampleMeshLandmarksFileList = new File(generalPath, "/Raw/fullMeshes/landmarks/noisyMeshLandmarks/").listFiles().toSeq

      (modelFile, modelLmsFile, sampleMeshFileList, sampleMeshLandmarksFileList)
  }
  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUI()
    val sampleFiles = new File(generalPath, "Registered/partialMeshes/ICPProposal").listFiles().toSeq
    val modelLmsFile = new File(generalPath, "Registered/model/AugmentedPCA.json")
    val targetMeshFile = new File(generalPath, "Raw/partialMeshes/mirroredMeshes/").listFiles().toSeq
    val targetLmsFile = new File(generalPath, "Raw/partialMeshes/landmarks/").listFiles().toSeq

    val modelLms = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get

    sampleFiles.map { s =>
      val sampleMesh = MeshIO.readMesh(s).get
      ui.show(sampleMesh,s"${s.getName}")

      val sampleMeshesLandmarks = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get
     //  ui.show(sampleMeshesLandmarks,s"landmarks")
      val targetMeshFile = new File(generalPath, "Raw/partialMeshes/mirroredMeshes/").listFiles().toSeq
      val targetLmsFile = new File(generalPath, "Raw/partialMeshes/landmarks/").listFiles().toSeq
      val name= "00"+s.getName.replace(".vtk",".ply")
      println(name)
      val target= targetMeshFile.find{t => t.getName() == name}.get
      val targetMesh2 = MeshIO.readMesh(target).get
      ui.show(targetMesh2,s"${name}")
      val nameLM= "00"+s.getName.replace(".vtk",".json")
      val targetLms= targetLmsFile.find{l => l.getName() == nameLM}.get
      val lm = LandmarkIO.readLandmarksJson[_3D](targetLms).get

      println(name)
      val alignTransform = AlignmentTransforms.computeTransform(lm, sampleMeshesLandmarks, Point3D(0, 0, 0))

      val alignedMesh = targetMesh2.transform(alignTransform)
      val alignedLandmarks = lm.map(_.transform(alignTransform))
      ui.show(alignedMesh,s"${s.getName}")
    }
  }
  }

