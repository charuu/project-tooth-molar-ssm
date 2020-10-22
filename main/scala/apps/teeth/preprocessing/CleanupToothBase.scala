package apps.teeth.preprocessing

import java.io.File

import scalismo.common.PointId
import scalismo.geometry.{Point, _3D}
import scalismo.io.MeshIO
import scalismo.mesh.{TriangleCell, TriangleMesh3D}
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}


object CleanupToothBase {

  def boundaryIds(m: TriangleMesh3D): IndexedSeq[PointId] = {
    m.pointSet.pointIds.filter(id => m.operations.pointIsOnBoundary(id)).toIndexedSeq
  }

  def ids2points(m: TriangleMesh3D, ids: IndexedSeq[PointId]): IndexedSeq[Point[_3D]] = {
    ids.map(id => m.pointSet.point(id))
  }

  def pointIdTriangles(m: TriangleMesh3D, pId: PointId): Seq[TriangleCell] = {
    m.triangulation.triangles.filter(f => f.ptId1 == pId || f.ptId2 == pId || f.ptId3 == pId)
  }

  def triangleCos(m: TriangleMesh3D, t: TriangleCell): Double = {
    val p1 = m.pointSet.point(t.ptId1)
    val p2 = m.pointSet.point(t.ptId2)
    val p3 = m.pointSet.point(t.ptId3)
    val a = (p1-p2).normalize
    val b = (p1-p3).normalize
    val c = (p2-p3).normalize
    val angles = Seq((a dot b), (a dot c), (b dot c)).map(math.abs)
    angles.min
  }

  def idTriangleMap(m: TriangleMesh3D): Map[PointId, Seq[Double]] = {
    m.pointSet.pointIds.toIndexedSeq.map{ id =>
      val triangles = pointIdTriangles(m, id)
      val areas = triangles.map(t => m.computeTriangleArea(t))
      (id, areas)
    }.toMap
  }

  def triangleAreas2ratio(areas: Seq[Double]): Double ={
    val avg = areas.sum / areas.length.toDouble
    (areas.map(a => a/avg).sum)/(areas.length).toDouble
  }

  def triangleAreas2minMax(areas: Seq[Double]): (Double, Double) ={
    val avg = areas.sum / areas.length.toDouble
    val ratios = areas.map(a => a/avg)
    (ratios.min, ratios.max)
  }


  @scala.annotation.tailrec
  def recursionPointIdAdd(target: TriangleMesh3D, allIds: Seq[PointId], newIds: Seq[PointId], maxIter: Int = 100): Seq[PointId] = {
    println(s"Iteration: ${maxIter}, allIds: ${allIds.length}")
    val newIdsfiltered = newIds.filter{id =>
      val triangles = pointIdTriangles(target, id)
      val areas = triangles.map(t => target.computeTriangleArea(t))
      val (ratioMin, ratioMax) = triangleAreas2minMax(areas)
      //      val angles = triangles.map{t => triangleCos(target, t)}
      val decision = ratioMin > 0.7 && ratioMax < 1.3 //&& angles.map(a => a < 0.7).foldLeft(true)(_ && _)
      decision
    }

    val newAllIds = (allIds ++ newIdsfiltered).toSet.toIndexedSeq
    val newnewIds = newIdsfiltered.flatMap(id => pointIdTriangles(target, id).flatMap(t => Seq(t.ptId1, t.ptId2, t.ptId3))).filter(id => !newAllIds.contains(id)).toSet.toIndexedSeq

    if(newnewIds.length < 5 || maxIter == 0) allIds
    else recursionPointIdAdd(target, newAllIds, newnewIds, maxIter-1)
  }

  def main(args: Array[String]): Unit = {
    println("Tooth cleaner!")
    scalismo.initialize()

    val inputPath = new File("/export/skulls/projects/teeth/data/surface/ok6er_extended")
    val outputPath = new File("/export/skulls/projects/teeth/data/surface/ok6er/")
    val meshes: Array[File] = inputPath.listFiles(f => f.getName.endsWith(".stl")).sorted
    //    val targetFile = meshes.slice(10, 11).head
    meshes.par.foreach { targetFile =>
      val target = MeshIO.readMesh(targetFile).get
      println(s"Target: ${targetFile}, Total number of points: ${target.pointSet.numberOfPoints}")
      val bIds = boundaryIds(target)
      val pointToRemove = recursionPointIdAdd(target, Seq(), bIds)

      val masked = target.operations.maskPoints(f => !pointToRemove.contains(f))
      val clipped = masked.transformedMesh
      val outputFile = new File(outputPath, targetFile.getName)
      MeshIO.writeMesh(clipped, outputFile)
      println(s"Clipped: ${targetFile}, Total number of points: ${clipped.pointSet.numberOfPoints}, output: ${outputFile}")

      val ui = ScalismoUIHeadless()
      ui.show(target, "target")
      ui.show(clipped, "masked")
    }
  }
}
