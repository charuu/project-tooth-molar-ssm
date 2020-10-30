package apps.teeth.utilities

import java.io.File

import Paths.generalPath
import apps.teeth.registeration.IcpRegistration.meshType
import apps.util.AlignmentTransforms
import breeze.linalg.{DenseMatrix, DenseVector, trace}
import scalismo.common.{DiscreteField, Domain, NearestNeighborInterpolator, PointId, ScalarArray, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector, Landmark, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatismoIO}
import scalismo.mesh.{MeshMetrics, ScalarMeshField, TriangleMesh}
import scalismo.numerics.Sampler
import scalismo.registration.GaussianProcessTransformationSpace
import scalismo.statisticalmodel.{MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

import scala.collection.immutable.Stream.Empty



case class MyLMSampler3D(mesh: TriangleMesh[_3D], lmPoints: Seq[Point[_3D]], numOfNeightbours: Int, numberOfPoints: Int)(implicit rng: Random)
  extends Sampler[_3D] {

  override val volumeOfSampleRegion = mesh.area

  private val p: Double = 1.0 / mesh.area

  val samplePointsInit: Set[Point[_3D]] = lmPoints.flatMap(lm => (mesh.pointSet.findNClosestPoints(lm, numOfNeightbours).map(_.point))).toSet
  val samplePoints: IndexedSeq[(Point[_3D], Double)] = scala.util.Random.shuffle(samplePointsInit).take(numberOfPoints).map(f => (f, p)).toIndexedSeq



  override def sample() = samplePoints
}

object  MyLMSampler3D {
  val ui = ScalismoUI()


 // val modelLmsFile2 = new File(generalPath, "/Ref/landmarks/ref_landmarks.json")
  //val modelLms2 = modelLms

 // val meshFile = new File(generalPath, "Registered/partialMeshes/0012_36good.vtk")
 // val mesh = MeshIO.readMesh(meshFile).get

 // val modelFile = new File(generalPath, "Registered/model/augmentedPcaModel/AugmentedPCA_0606.h5")
 // val model = StatismoIO.readStatismoMeshModel(modelFile).get
 // val modelLmsFile = new File(generalPath, "Registered/model/landmarks/AugmentedPCA.json")
 // val modelLms = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get

 // val targetMeshFile = new File(generalPath, "Raw/partialMeshes/mirroredMeshes/0012_36.ply")
 // val targetMeshInit = MeshIO.readMesh(targetMeshFile).get
//  val targetLmsFile = new File(generalPath, s"Raw/partialMeshes/landmarks/0012_36.json")

  //val targetLmsInit = LandmarkIO.readLandmarksJson[_3D](targetLmsFile).get

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    implicit val rng = scalismo.utils.Random(42)
    val (model, modelLms, targetMeshes, targetLms) = LoadTestData.modelAndTarget(meshType)
    val alignedMeshLms = alignShapes(modelLms,targetMeshes,targetLms).align()
    Viewer(model,model.referenceMesh,targetMeshes.seq(0)._1,"sampledTarget")



  }

  case class Viewer(model: StatisticalMeshModel,mesh:TriangleMesh[_3D],fixedMesh:TriangleMesh[_3D],name:String){

    val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
    val lowRankGP = model.gp.interpolate(interpolator)
    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)
    val observationGroup = ui.createGroup("observation")


    val seqModelCorr = fixedMesh.pointSet.points.toIndexedSeq.map{
      pt=> (pt,mesh.pointSet.findClosestPoint(mesh.operations.closestPointOnSurface(pt).point).id)
    }

    val projection = (pt : Point[_3D]) => {
    val fixedPo = seqModelCorr.find( p => mesh.pointSet.findClosestPoint(pt).id.id == p._2.id)

      if(fixedPo != None) {
       fixedMesh.operations.closestPointOnSurface(fixedPo.get._1).point
      } else
        mesh.pointSet.findClosestPoint(pt).point}


    val registrationTransformation = transformationSpace.transformForParameters(model.coefficients(mesh))
    val finalTransformation = registrationTransformation.andThen( projection)

    val projectedMesh = mesh.transform(finalTransformation)

    val resultGroup = ui.createGroup("result")
    val projectionView = ui.show(resultGroup, projectedMesh, "projection")


  }
}


