package apps.teeth.utilities

import java.io.File

import Paths.generalPath
import apps.util.AlignmentTransforms
import breeze.linalg.{DenseMatrix, DenseVector, trace}
import scalismo.common.{DiscreteField, NearestNeighborInterpolator, PointId, ScalarArray, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector, Landmark, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatismoIO}
import scalismo.mesh.{MeshMetrics, ScalarMeshField, TriangleMesh}
import scalismo.numerics.Sampler
import scalismo.registration.GaussianProcessTransformationSpace
import scalismo.statisticalmodel.{MultivariateNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random



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
  val modelLmsFile2 = new File(generalPath, "/Ref/landmarks/ref_landmarks_7.json")
  val modelLms2 = LandmarkIO.readLandmarksJson[_3D](modelLmsFile2).get

  val meshFile = new File(generalPath, "registeredPartialMeshes/Target.vtk")
  val mesh = MeshIO.readMesh(meshFile).get

  val modelFile = new File(generalPath, "Registered/model/AugmentedPCA_0606.h5")
  val model = StatismoIO.readStatismoMeshModel(modelFile).get
  val modelLmsFile = new File(generalPath, "Registered/model/AugmentedPCA.json")
  val modelLms = LandmarkIO.readLandmarksJson[_3D](modelLmsFile).get

  val targetMeshFile = new File(generalPath, "Raw/partial/mirrored_meshes/0021_36.ply")
  val targetMeshInit = MeshIO.readMesh(targetMeshFile).get
  val targetLmsFile = new File(generalPath, s"Raw/partial/landmarks/0021_36.json")

  val targetLmsInit = LandmarkIO.readLandmarksJson[_3D](targetLmsFile).get

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    implicit val rng = scalismo.utils.Random(42)


    val lmNames = Seq("F1", "F2","L1","L2")
    val targetTransform = AlignmentTransforms.computeTransform( targetLmsInit, modelLms,Point(0,0,0))
    val lmsTransform = targetLmsInit.map(_.transform(targetTransform))

    val lmPoints = lmsTransform.filter(lm => lmNames.contains(lm.id)).map(_.point)


    val fixedMesh = targetMeshInit.transform(targetTransform)
    ui.show(fixedMesh,"fixedMesh")
    val mySampler = MyLMSampler3D(fixedMesh, lmPoints, 50000, 15000)

    Viewer(model,mesh,fixedMesh,mySampler,"sampledTarget")



  }

  case class Viewer(model: StatisticalMeshModel,mesh:TriangleMesh[_3D],fixedMesh:TriangleMesh[_3D], mySampler: MyLMSampler3D,name:String){

    val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
    val lowRankGP = model.gp.interpolate(interpolator)
    val transformationSpace = GaussianProcessTransformationSpace(lowRankGP)
    val observationGroup = ui.createGroup("observation")


    val projection = (pt : Point[_3D]) => {

      val meshid =  mesh.operations.closestPointOnSurface(pt)
      val m = mesh.pointSet.findClosestPoint(meshid.point)
      val fixedid= fixedMesh.pointSet.findClosestPoint(pt)
      val avgDist = MeshMetrics.avgDistance(fixedMesh,mesh)
   //   println(math.(avdDist,)avdDist)
      val bool = Math.round(fixedMesh.operations.shortestDistanceToSurfaceSquared(meshid.point))
      val intersection = fixedMesh.operations.getIntersectionPointsOnSurface(m.point,mesh.vertexNormals.atPoint(m.id)).length
    //  println(BigDecimal(Math.sqrt(bool)).setScale(2,BigDecimal.RoundingMode.HALF_UP) )
       if( intersection !=0 && bool <= 1) {
         fixedid.point
         }

      else{
         m.point
        }
    }
    val registrationTransformation = transformationSpace.transformForParameters(model.coefficients(mesh))
    val finalTransformation = registrationTransformation.andThen( projection)

    val projectedMesh = model.mean.transform(finalTransformation)

    val resultGroup = ui.createGroup("result")
    val projectionView = ui.show(resultGroup, projectedMesh, "projection")


  //  val uncert2:ScalarArray[Double] = ScalarArray[Double](posteriorFromTarget.mean.pointSet.pointIds.toIndexedSeq.map{id=> trace(posteriorFromTarget.gp.cov(id,id)) }.toArray)
  //  val cMesh2: ScalarMeshField[Double] = ScalarMeshField(posteriorFromTarget.mean,uncert2)

  //  ui.show(cMesh2,"Mesh ")
  }
}


