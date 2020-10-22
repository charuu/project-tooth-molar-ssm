package apps.teeth.utilities

import java.io.File

import Paths.generalPath
import apps.util.{AlignmentTransforms, FileUtils}
import scalismo.geometry.{Landmark, Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.utils.Random

object AlignShapes {

  def main(args: Array[String]) {
    scalismo.initialize()
    val (model, modelLms, sampleMeshes, sampleLMs) = LoadTestData.modelAndTarget("Partial")

    alignShapes(modelLms,sampleMeshes,sampleLMs)

  }
}

case class alignShapes(modelLms:Seq[Landmark[_3D]], sampleMeshes:Seq[(TriangleMesh[_3D],String)], sampleLms:Seq[(Seq[Landmark[_3D]],String)]) {
  implicit val random: Random = Random(1024)
  private val origin = Point3D(0, 0, 0)
  private val alignedMeshesPath = new File(generalPath, "/aligned/meshes")
  private val alignedLMsPath = new File(generalPath, "/aligned/landmarks")

  def align(): Seq[((TriangleMesh[_3D],String),(Seq[Landmark[_3D]],String))] ={
    sampleMeshes.map { s =>
      val (mesh , meshName) = (s._1 ,s._2)
      val (landmarks , lmName) = sampleLms.find(l => l._2 == meshName.replace(".ply",".json")).get

      val alignTransform = AlignmentTransforms.computeTransform(landmarks, modelLms, origin)

      val alignedMesh = mesh.transform(alignTransform)
      val alignedLandmarks = landmarks.map(_.transform(alignTransform))

      MeshIO.writeVTK(alignedMesh, new File(alignedMeshesPath, meshName))
      LandmarkIO.writeLandmarksJson[_3D](alignedLandmarks, new File(alignedLMsPath, lmName))

      ((alignedMesh,meshName),(alignedLandmarks,lmName))
    }
  }
}