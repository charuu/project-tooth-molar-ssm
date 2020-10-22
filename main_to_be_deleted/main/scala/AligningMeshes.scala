/*
 *  Copyright University of Basel, Graphics and Vision Research Group
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

import java.io.File

import scalismo.geometry.{Landmark, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.registration.LandmarkRegistration
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.ScalismoUI

import main.scala.Paths.generalPath
object AligningMeshes {
 scalismo.initialize()
 implicit val rng = scalismo.utils.Random(42)
 val ui = ScalismoUI()
 val model = StatisticalModelIO.readStatisticalMeshModel(new java.io.File(generalPath + "/Ref/model/tooth_gp_model_1000-components.h5")).get
///Registered/model/AugmentedPCA.h5"

 def main(args: Array[String]): Unit = {


  val modelLm = LandmarkIO.readLandmarksJson[_3D](new java.io.File(generalPath +"/Ref/landmarks/ref_landmarks.json")).get
  val referenceLandmarkViews = modelLm.map(lm => ui.show( lm, s"lm-${lm.id}"))
  ui.show(model,"Model")
  val sampleLmFileList = new java.io.File(generalPath +"/Raw/landmarks/").listFiles.sortBy(f => f.getName())

  val sampleMeshesFileList = new java.io.File(generalPath +"/Raw/mesh").listFiles.sortBy(f => f.getName())
  val experimentRegistrationPath = new File(generalPath +"/aligned/meshes/")
  val experimentRegistrationPath2 = new File(generalPath +"/aligned/landmarks/")
  val targetGroup = ui.createGroup("target")

  sampleMeshesFileList.foreach { f =>
   val target = MeshIO.readMesh(f).get
   val flag = new java.io.File(generalPath + "/Raw/landmarks/"+ f.getName.replace(".stl", ".json")).exists
   if (flag == true) {

    val targetLM = LandmarkIO.readLandmarksJson[_3D](new java.io.File(generalPath + "/Raw/landmarks/"+ f.getName.replace(".stl", ".json"))).get
    targetLM.map(lm => ui.show(targetGroup, lm, s"lm-${lm.id}"))

    val alignedTarget = alignModel(model, target, targetLM, modelLm)
    ui.show(alignedTarget._1, "aligned")
    MeshIO.writeMesh(alignedTarget._1, new File(experimentRegistrationPath, f.getName))
    LandmarkIO.writeLandmarksJson(alignedTarget._2, new File(experimentRegistrationPath2, f.getName.replace(".stl", ".json")))
   }
  }
 }

 def alignModel(model: StatisticalMeshModel, croppedMesh: TriangleMesh[_3D], targetLandmark: Seq[Landmark[_3D]], modelLandmark: Seq[Landmark[_3D]]): (TriangleMesh[_3D],Seq[Landmark[_3D]]) = {
  val bestTransform = LandmarkRegistration.rigid3DLandmarkRegistration( modelLandmark,targetLandmark, center = Point(0, 0, 0))
  val aligned = model.transform(bestTransform)
  val transLM = modelLandmark.map { lm => lm.copy(point = bestTransform(lm.point)) }.sortBy(lm => lm.id)

  (aligned.mean,transLM)
 }
}