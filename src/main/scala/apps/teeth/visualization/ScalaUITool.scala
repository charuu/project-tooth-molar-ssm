package apps.teeth.visualization

import java.awt.Color
import java.io.File

import api.sampling.proposals.MyLMSampler3D
import scalismo.common.{NearestNeighborInterpolator, PointId, UnstructuredPointsDomain}
import scalismo.utils.Random

//import apps.teeth.registeration.NonRigidRegistration.ui
import apps.teeth.utilities.Paths.generalPath
import apps.teeth.utilities.{LoadTestData, alignShapes}
import apps.util.{AlignmentTransforms, FileUtils}
import breeze.linalg.trace
import scalismo.common.ScalarArray
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{ImageIO, LandmarkIO, MeshIO, StatismoIO}
import scalismo.mesh.{MeshMetrics, ScalarMeshField}
import scalismo.ui.api.ScalismoUI

object ScalaUITool {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUI()
    implicit val random: Random = Random(1024)
    val (m, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget("Full")
 //   val alignedMeshLms = alignShapes(modelLms,targetMeshes,targetLms).align()
    val modelGroup = ui.createGroup("model")
    val model = m.decimate((5000))
    ui.show(modelGroup,model,"model")
    val pcaLandmark = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "/Ref/landmarks/ref_landmarks.json")).get
 ui.show(pcaLandmark(0),"lm").color = Color.RED
   val modelPoints = MyLMSampler3D(model.referenceMesh,Seq(pcaLandmark(0).point), 1300,1300).sample().map(_._1)
     val modelIds: IndexedSeq[PointId] = modelPoints.map(p => model.referenceMesh.pointSet.findClosestPoint(p).id).toIndexedSeq

     val showModel = ui.show(modelGroup, UnstructuredPointsDomain(modelPoints), "model").radius = 0.1

   /* alignedMeshLms.map { targetMeshLms =>


        val target = targetMeshLms._1._1
        val targetName = targetMeshLms._1._2
        val targetLms = targetMeshLms._2._1

        val targetGroup = ui.createGroup("target")
        ui.show(targetGroup, target, s"$targetName")
      ui.show(targetGroup,targetLms,"lm")

      } */
/*

    ui.createGroup("Mesh")
  //  val modelFile = new File(generalPath, "/Ref/model/tooth_gp_model_50000-components.h5")
    val modelFile = new File(generalPath, "Registered/model/pcaModel_1000points_0406.h5")
    val m = StatismoIO.readStatismoMeshModel(modelFile).get

    val uncert:ScalarArray[Double] = ScalarArray[Double](m.referenceMesh.pointSet.pointIds.toIndexedSeq.map{id=> trace(m.gp.cov(id,id)) }.toArray)
    val cMesh: ScalarMeshField[Double] = ScalarMeshField(m.referenceMesh,uncert)

    ui.show(cMesh,"pca ")
    val modelFile2 = new File(generalPath, "Registered/model/AugmentedPCA_0406B.h5")
    val m2 = StatismoIO.readStatismoMeshModel(modelFile2).get

    val uncert2:ScalarArray[Double] = ScalarArray[Double](m2.referenceMesh.pointSet.pointIds.toIndexedSeq.map{id=> trace(m2.gp.cov(id,id)) }.toArray)
    val cMesh2: ScalarMeshField[Double] = ScalarMeshField(m2.referenceMesh,uncert2)

    ui.show(cMesh2,"augmented pca ")

    val sampleMeshesFileList = new java.io.File(generalPath+"/10LowerMolarData/rawMeshes/").listFiles.sortBy(f => f.getName())
    val croppedimagesFileList = new java.io.File(generalPath+"/10LowerMolarData/croppedImages/").listFiles.sortBy(f => f.getName())
    val sampleMeshesRegisteredFileList = new java.io.File(generalPath+"/10LowerMolarData/registered/").listFiles.sortBy(f => f.getName())

    val initialPath = new File(generalPath, "/Raw")
    val alignedPath = new File(generalPath, "/aligned")
    val initialLMsPath = new java.io.File(initialPath + "/landmarks/noisyMeshLandmarks/")
    val alignedLMsPath = new File(alignedPath + "/landmarks")
    val origin = Point3D(0, 0, 0)

    sampleMeshesFileList.map{
      m =>
        val basename = FileUtils.basename(m)
        val lmName = s"$basename.json"
        val mesh = MeshIO.readMesh(m).get
        val flag = new File(initialLMsPath, lmName).exists
        if(flag ==true) {


          val grp = ui.createGroup(s"${m.getName}")
        //  val mAligned = MeshIO.readMesh(m).get

             val mReg =   MeshIO.readMesh(new java.io.File(generalPath+"/10LowerMolarData/registered/" + m.getName.replace(".stl",".vtk"))).get
          val lms = LandmarkIO.readLandmarksJson[_3D](new File(initialLMsPath, lmName)).get
          val alms = LandmarkIO.readLandmarksJson[_3D](new File(alignedLMsPath, lmName)).get
          val transform = AlignmentTransforms.computeTransform(alms, lms, origin)
          val tm = mReg.transform(transform)
          ui.show(grp, mesh, s"Mesh ${m.getName}")

           ui.show(grp,tm,s"Registered- ${m.getName} -${MeshMetrics.avgDistance(mesh,tm)}")
          val image = ImageIO.read3DScalarImage[Short](new java.io.File(generalPath + "/10LowerMolarData/croppedImages/" + m.getName.replace(".stl", "_cropped.nii"))).get
          ui.show(grp, image, s"image- ${m.getName}")
        }
    }
*/

  }
}
