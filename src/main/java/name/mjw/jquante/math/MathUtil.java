package name.mjw.jquante.math;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import name.mjw.jquante.math.geom.Point3D;

/**
 * A collection of few misc. utility math functions. All methods are static and
 * the class cannot be instantiated.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class MathUtil {

	/** Creates a new instance of MathUtil */
	private MathUtil() {
	}

	/**
	 * Method to convert radians to degrees
	 * 
	 * @param radians
	 *            - the value
	 * @return the equivalent in degrees
	 */
	public static double toDegrees(double radians) {
		return (radians * 180.0 / Math.PI);
	}

	/**
	 * Method to convert degrees to radians
	 * 
	 * @param degrees
	 *            - the value
	 * @return the equivalent in radians
	 */
	public static double toRadians(double degrees) {
		return (degrees * Math.PI / 180.0);
	}

	/**
	 * @param val
	 *            a double value
	 * @return true if val is within 2e-6 of zero
	 */
	public static boolean isNearZero(double val) {
		return isNearZero(val, 2e-6);
	}

	/**
	 * 
	 * @param val
	 *            a double value
	 * @param epsilon
	 *            the "near" distance from zero
	 * @return true if val is within a distance epsilon of zero
	 */
	public static boolean isNearZero(double val, double epsilon) {
		return (Math.abs(val) < epsilon);
	}

	/**
	 * Method to find the angle in radians defined by three points v1-v2-v3.
	 * 
	 * @param v1
	 *            the first point
	 * @param v2
	 *            the second point (central angle)
	 * @param v3
	 *            the third point
	 * @return the angle defined
	 */
	public static double findAngle(Point3D v1, Point3D v2, Point3D v3) {
		Vector3D v12 = new Vector3D(v2.sub(v1));
		Vector3D v32 = new Vector3D(v2.sub(v3));

		return v12.angleWith(v32);
	}

	/**
	 * Method to find the dihedral angle defined by planes v1-v2-v3 and
	 * v2-v3-v4.
	 * 
	 * @param v1
	 *            first point
	 * @param v2
	 *            second point
	 * @param v3
	 *            third point (2nd and 3rd point define the angle)
	 * @param v4
	 *            the fourth angle
	 * @return the dihedral angle defined
	 */
	public static double findDihedral(Point3D v1, Point3D v2, Point3D v3,
			Point3D v4) {
		// normal of plane 1
		Vector3D v12 = new Vector3D(v2.sub(v1));
		Vector3D v32 = new Vector3D(v2.sub(v3));
		Vector3D n123 = v12.cross(v32).normalize();

		// normal of plane 2
		Vector3D v23 = new Vector3D(v3.sub(v2));
		Vector3D v43 = new Vector3D(v3.sub(v4));
		Vector3D n234 = v23.cross(v43).normalize();

		// sign of the dihedral
		double sign = v32.mixedProduct(n123, n234);

		if (sign >= 0.0)
			sign = -1.0;
		else
			sign = 1.0;

		// and find the angle between the two planes
		return n123.angleWith(n234) * sign;
	}

	/**
	 * compute N!
	 * 
	 * @param n
	 *            the n, whose factorial is to be found
	 * @return the factorial
	 */
	public static int factorial(int n) {
		int value = 1;

		while (n > 1) {
			value = value * n;
			n--;
		}
		return value;
	}

	/**
	 * compute double N! ... (1*3*5*...*n)
	 * 
	 * @param n
	 *            the n, whose factorial is to be found
	 * @return the factorial
	 */
	public static int factorial2(int n) {
		int value = 1;

		while (n > 0) {
			value = value * n;
			n -= 2;
		}
		return value;
	}

	/**
	 * Does ( a! / b! / (a-2*b)! )
	 * 
	 * @param a
	 *            the first term
	 * @param b
	 *            the second term
	 * @return ( a! / b! / (a-2*b)! )
	 */
	public static double factorialRatioSquared(int a, int b) {
		return factorial(a) / factorial(b) / factorial(a - 2 * b);
	}

	/**
	 * Pre-factor of binomial expansion.
	 * 
	 * From Augspurger and Dykstra: <a
	 * href="http://dx.doi.org/10.1021/j100176a037">10.1021/j100176a037</a>
	 * 
	 * @param s
	 *            s
	 * @param ia
	 *            ia
	 * @param ib
	 *            ib
	 * @param xpa
	 *            xpa
	 * @param xpb
	 *            xpb
	 * @return Pre-factor of binomial expansion.
	 */
	public static double binomialPrefactor(int s, int ia, int ib, double xpa,
			double xpb) {
		double sum = 0.0;

		for (int t = 0; t < (s + 1); t++) {
			if (((s - ia) <= t) && (t <= ib)) {
				sum += binomial(ia, s - t) * binomial(ib, t)
						* Math.pow(xpa, ia - s + t) * Math.pow(xpb, ib - t);
			}
		}

		return sum;
	}

	/**
	 * Binomial coefficient
	 * 
	 * Number of ways to choose k elements from a set of n elements.
	 * 
	 * @param n
	 *            number of n elements in set
	 * @param k
	 *            number of k elements in set
	 * @return binomial coefficient
	 * 
	 *         https://en.wikipedia.org/wiki/Binomial_coefficient
	 */
	public static int binomial(int n, int k) {
		return (factorial(n) / factorial(k) / factorial(n - k));
	}

	/**
	 * Relative difference between two floating point numbers, as described at:
	 * <a href="http://c-faq.com/fp/fpequal.html">
	 * http://c-faq.com/fp/fpequal.html</a>
	 * 
	 * @param a
	 *            first number
	 * @param b
	 *            second number
	 * @return the relative difference
	 */
	public static double relativeDifference(double a, double b) {
		double c = Math.abs(a);
		double d = Math.abs(b);

		d = Math.max(c, d);

		return d == 0.0 ? 0.0 : Math.abs(a - b) / d;
	}

	public static RealVector RealMatrixToRealVector(RealMatrix realMatrix) {
		RealVector vec = new ArrayRealVector(realMatrix.getRowDimension() * realMatrix.getColumnDimension());
		int ii = 0;

		for (int i = 0; i < realMatrix.getRowDimension(); i++)
			for (int j = 0; j < realMatrix.getColumnDimension(); j++)
				vec.setEntry(ii++, realMatrix.getEntry(i, j));

		return vec;

	}
}