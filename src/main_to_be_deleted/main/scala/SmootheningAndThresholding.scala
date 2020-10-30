import java.io.File

import main.scala.Paths.generalPath
import scalismo.geometry._3D
import scalismo.image.DiscreteScalarImage
import scalismo.io.{ImageIO, MeshIO}
import scalismo.mesh.{TriangleMesh, TriangleMesh3D}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.{ImageConversion, MeshConversion}
import vtk._

object SmootheningAndThresholding {
  val ui = ScalismoUI()
  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val file = new File(generalPath +"/Raw/croppedimages/2.5.90_cropped.nii")
    val img: DiscreteScalarImage[_3D, Short] = ImageIO.read3DScalarImageAsType[Short](file, resampleOblique = true).get
    val mesh = extractContour(img, 800) // 250 usually
    MeshIO.writeMesh(mesh, new File(generalPath+"/Raw/croppedimages/2.5.90.stl"))
    ui.show(mesh,"mesh")
  }

  def getBiggestConnectedComponent(mesh: TriangleMesh3D, componentThreshold: Int): TriangleMesh3D = {
    val connectivity = new vtkPolyDataConnectivityFilter()

    val polydata = MeshConversion.meshToVtkPolyData(mesh)

    connectivity.SetExtractionModeToSpecifiedRegions()
    connectivity.SetInputData(polydata)
    connectivity.Update()
    val count = connectivity.GetNumberOfExtractedRegions()

    val idListFull = {
      val vtk = connectivity.GetRegionSizes()
      val out = for (i <- -0 until count) yield vtk.GetValue(i)
      vtk.Delete()
      out
      }.zipWithIndex.sortBy(_._1).reverse

    connectivity.InitializeSpecifiedRegionList()

    val idList = idListFull.map(id => id._2)

    for(i <- idListFull.indices if idListFull(i)._1 > componentThreshold){
      connectivity.AddSpecifiedRegion(idListFull(i)._2)
    }

    connectivity.Update()

   /* val fill = new vtkFillHolesFilter()
    fill.SetInputConnection(connectivity.GetOutputPort())
    fill.SetHoleSize(1.0)
    fill.Update()*/

    val out = new vtkPolyData()
    out.DeepCopy(connectivity.GetOutput())
    MeshConversion.vtkPolyDataToCorrectedTriangleMesh(out).get

  }

  def extractContour(img : DiscreteScalarImage[_3D, Short], thresholdValueLB: Double = 0.5): TriangleMesh3D = {
    //    DiscreteScalarImage()
    println("Thresholds set to " + thresholdValueLB)
    val im = img.map(value => if (value > thresholdValueLB) 1 else 0)
    ui.show(im,"image")
    println("turning into vtkStructuredPoints...")
    val imgvtk: vtkStructuredPoints = ImageConversion.imageToVtkStructuredPoints(im)

    println("making new vtkContourFilter...")
    val mc = new vtkContourFilter()
    mc.SetValue(0, 0.01)
    mc.ComputeNormalsOn()
    mc.ComputeGradientsOn()
    mc.SetNumberOfContours(1)
    mc.SetInputData(imgvtk)


    println("decimating...")
    val reduction = 0.0 // Reduction value sets the % that should be "removed" 0.2 means reduce to 80% of original
    val decimate = new vtkDecimatePro()
    println("new decimate done")
    decimate.SetInputConnection(mc.GetOutputPort())
    println("set input connection done")
    decimate.SetTargetReduction(reduction)
    println("set target reduction done")
    decimate.Update()
    println("update decimate done")

    println("smoothing...")
    val smooth = new vtkSmoothPolyDataFilter
    smooth.SetInputConnection(decimate.GetOutputPort())
    smooth.SetNumberOfIterations(50)
    smooth.SetRelaxationFactor(0.1)
    smooth.FeatureEdgeSmoothingOn()
    smooth.SetFeatureAngle(120)
    smooth.BoundarySmoothingOn()
    smooth.Update()
    smooth.GetConvergence()

    println("finding normals...")
    val normalGenerator = new vtkPolyDataNormals
    normalGenerator.SetInputConnection(smooth.GetOutputPort())
    normalGenerator.ComputePointNormalsOn()
    normalGenerator.ComputeCellNormalsOn()
    normalGenerator.Update()

    println("getting vtkMesh...")
    val meshVTK: vtkPolyData = normalGenerator.GetOutput()

    mc.Delete()
    decimate.Delete()
    smooth.Delete()
    normalGenerator.Delete()

    println("getting triangle mesh...")
    val surface: TriangleMesh[_3D] = MeshConversion.vtkPolyDataToCorrectedTriangleMesh(meshVTK).get
    ui.show(surface,"show")
    println("deleting meshVTK...")
    meshVTK.Delete()

    println("saving biggest connected component...")
    // the image labeling is a bit messy leading to some floating points/triangles outside of the mesh
    // as these can mess up bot the registration and the mesh metrics, we extract the largest regions
    val s = getBiggestConnectedComponent(surface, 7000)

   // MeshIO.writeMesh(s, new File("data/LowerMolar/Raw/croppedimages/mesh.stl"))
   // ui.show(s,"show")
    s
  }
}
