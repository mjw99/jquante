package name.mjw.jquante;

import javax.vecmath.Point3d;

public class GaussianTypeOrbital {

	/*
	 * g(x,y,z) = A*(x^i)*(y^j)*(z^k)*exp{-a*(r-ro)^2}
	 */

	private double exponent;

	private Point3d origin = new Point3d(0, 0, 0);

	/**
	 * Quantum number l
	 */
	private int l = 0;
	/**
	 * Quantum number m
	 */
	private int m = 0;
	/**
	 * Quantum number n
	 */
	private int n = 0;

	private double norm = 1.0;

	public GaussianTypeOrbital(double exponent) {

		this.exponent = exponent;

		normalise();

	}

	public GaussianTypeOrbital(double exponent, Point3d origin, int l, int m,
			int n) {

		this.exponent = exponent;

		this.origin = origin;

		this.l = l;
		this.m = m;
		this.n = n;

		normalise();

	}

	private void normalise() {

		norm = Math.sqrt(Math.pow(2, 2 * (l + m + n) + 1.5)
				* Math.pow(exponent, l + m + n + 1.5)
				/ Utils.doubleFactorial(2 * l - 1)
				/ Utils.doubleFactorial(2 * m - 1)
				/ Utils.doubleFactorial(2 * n - 1) / Math.pow(Math.PI, 1.5));

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

		return norm * Math.pow(dx, l) * Math.pow(dx, m) * Math.pow(dx, n)
				* Math.exp(-1 * exponent * d2);

	}
}
