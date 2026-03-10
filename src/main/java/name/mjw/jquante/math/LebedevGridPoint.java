package name.mjw.jquante.math;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;

import net.jafama.FastMath;

/**
 * Point information on a Lebedev grid.
 *
 * @param x      x coordinate of the point in the unit sphere
 * @param y      y coordinate of the point in the unit sphere
 * @param z      z coordinate of the point in the unit sphere
 * @param weight weight of the point
 */
public record LebedevGridPoint(double x, double y, double z, double weight) {

	/**
	 * Returns the grid point as a Vector3D.
	 *
	 * @return the grid point
	 */
	public Vector3D getVector3D() {
		return new Vector3D(x, y, z);
	}

	/**
	 * Returns the grid point as a RealVector.
	 *
	 * @return the grid point
	 */
	public RealVector getXYZ() {
		return new ArrayRealVector(new double[] { x, y, z }, false);
	}

	/**
	 * Returns the grid point coelevation.
	 *
	 * @return the grid point coelevation in [-PI,PI]
	 */
	public double getPhi() {
		return FastMath.acos(z);
	}

	/**
	 * Returns the grid point azimuth.
	 *
	 * @return the grid point azimuth in [-PI,PI]
	 */
	public double getTheta() {
		return FastMath.atan2(y, x);
	}
}
