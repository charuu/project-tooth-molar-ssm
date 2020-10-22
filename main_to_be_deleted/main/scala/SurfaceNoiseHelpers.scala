package main.scala

import breeze.linalg.DenseVector
import scalismo.geometry._
import scalismo.statisticalmodel.{MultivariateNormalDistribution}

object SurfaceNoiseHelpers {

  /**
   * Calculates the multivariate normal distribution for surface normal dependant noise.
   *
   * @param surfaceNormal                 Normal of the surface at the given position.
   * @param stdDevNoiseAlongSurfaceNormal Standard deviation of the noise along the normal.
   * @param stdDevNoiseInTangentialPlane  Standard deviation of the noise in the tangent plane.
   */
  def surfaceNormalDependantNoise(
                                   surfaceNormal: EuclideanVector3D,
                                   stdDevNoiseAlongSurfaceNormal: Double,
                                   stdDevNoiseInTangentialPlane: Double
                                 ): MultivariateNormalDistribution = {

    val normalDirection = surfaceNormal.normalize

    /* There is an infinite number of perpendicular vectors, any one will do.
    * To find any perpendicular vector, just take the cross product with any other, non-parallel vector.
    * We try (1,0,0), and if it happened to be parallel, the crossproduct is (0,0,0), and we take another.
    */
    val firstTangentialDirection = {
      val candidate = normalDirection.crossproduct(EuclideanVector3D(1, 0, 0))
      val direction = if (candidate.norm2 < 0.0001) candidate else normalDirection.crossproduct(EuclideanVector3D(0, 1, 0))
      direction.normalize
    }

    val secondTangentialDirection = normalDirection.crossproduct(firstTangentialDirection).normalize

    val normalNoiseVariance = stdDevNoiseAlongSurfaceNormal * stdDevNoiseAlongSurfaceNormal
    val tangentialNoiseVariance = stdDevNoiseInTangentialPlane * stdDevNoiseInTangentialPlane

    MultivariateNormalDistribution(DenseVector.zeros[Double](3), Seq(
      (normalDirection.toBreezeVector, normalNoiseVariance),
      (firstTangentialDirection.toBreezeVector, tangentialNoiseVariance),
      (secondTangentialDirection.toBreezeVector, tangentialNoiseVariance)
    ))
  }
}