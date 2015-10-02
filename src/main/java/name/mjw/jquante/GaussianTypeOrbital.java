package name.mjw.jquante;

import javax.vecmath.Point3d;

/**
 * Gaussian Type Orbitial based upon the following function:
 * 
 * g(x,y,z) = A*(x^i)*(y^j)*(z^k)*exp{-a*(r-ro)^2}
 */
public class GaussianTypeOrbital {

	private double exponent;

	private Point3d origin = new Point3d(0, 0, 0);

	private int i = 0;
	private int j = 0;
	private int k = 0;

	private double norm = 1.0;

	public GaussianTypeOrbital(double exponent) {

		this.exponent = exponent;
		normalise();

	}

	public GaussianTypeOrbital(double exponent, Point3d origin, int l, int m,
			int n) {

		this.exponent = exponent;
		this.origin = origin;

		this.i = l;
		this.j = m;
		this.k = n;

		normalise();

	}

	private void normalise() {

		norm = Math.sqrt(Math.pow(2, 2 * (i + j + k) + 1.5)
				* Math.pow(exponent, i + j + k + 1.5)
				/ Utils.doubleFactorial(2 * i - 1)
				/ Utils.doubleFactorial(2 * j - 1)
				/ Utils.doubleFactorial(2 * k - 1) / Math.pow(Math.PI, 1.5));

	}

	public void overlap(GaussianTypeOrbital Other) {

	}

	/**
	 * Compute the amplitude of the primitive Gaussian at point x,y,z
	 * 
	 * @param point
	 */
	public double valueAtPoint(Point3d point) {

		double dx = point.x - origin.x;
		double dy = point.y - origin.y;
		double dz = point.z - origin.z;

		double d2 = Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2);

		return norm * Math.pow(dx, i) * Math.pow(dx, j) * Math.pow(dx, k)
				* Math.exp(-1 * exponent * d2);

	}
}
