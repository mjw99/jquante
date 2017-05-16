package name.mjw.jquante.math.qm.integral;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.Vector3D;
import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.math.qm.basis.Power;

/**
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class NuclearTerm implements IntegralsPackage {

	/** Creates an instance of NuclearTerm */
	public NuclearTerm() {
	}

	/**
	 * The nuclear attraction term.
	 * 
	 * <i> Taken from THO eq. 2.15 <i>
	 * 
	 * @param a
	 *            the coefficient of primitive Gaussian a.
	 * @param norm1
	 *            the normalization factor of primitive Gaussian a.
	 * @param power1
	 *            the orbital powers of primitive Gaussian a.
	 * @param alpha1
	 *            the coefficient of primitive Gaussian a.
	 * @param b
	 *            the coefficient of primitive Gaussian b.
	 * @param norm2
	 *            the normalization factor of primitive Gaussian b.
	 * @param power2
	 *            the orbital powers of primitive Gaussian b.
	 * @param alpha2
	 *            the coefficient of primitive Gaussian b.
	 * @param c
	 *            the location of nuclear centre.
	 * @return the nuclear attraction integral.
	 */
	public double nuclearAttraction(Point3D a, double norm1, Power power1,
			double alpha1, Point3D b, double norm2, Power power2,
			double alpha2, Point3D c) {

		Point3D product = IntegralsUtil.gaussianProductCenter(alpha1, a,
				alpha2, b);

		double rABSquared = a.distanceSquaredFrom(b);
		double rCPSquared = c.distanceSquaredFrom(product);

		double gamma = alpha1 + alpha2;

		double[] ax = constructAArray(power1.getL(), power2.getL(),
				product.getX() - a.getX(), product.getX() - b.getX(),
				product.getX() - c.getX(), gamma);

		double[] ay = constructAArray(power1.getM(), power2.getM(),
				product.getY() - a.getY(), product.getY() - b.getY(),
				product.getY() - c.getY(), gamma);

		double[] az = constructAArray(power1.getN(), power2.getN(),
				product.getZ() - a.getZ(), product.getZ() - b.getZ(),
				product.getZ() - c.getZ(), gamma);

		double sum = 0.0;
		int i, j, k;
		for (i = 0; i < ax.length; i++) {
			for (j = 0; j < ay.length; j++) {
				for (k = 0; k < az.length; k++) {
					sum += ax[i]
							* ay[j]
							* az[k]
							* IntegralsUtil.computeFGamma(i + j + k, rCPSquared
									* gamma);
				}
			}
		}

		return (-norm1 * norm2 * 2.0 * Math.PI / gamma
				* Math.exp(-alpha1 * alpha2 * rABSquared / gamma) * sum);
	}

	/**
	 * <i> THO eq. 2.18 and 3.1 <i>
	 */
	public double[] constructAArray(int l1, int l2, double pa, double pb,
			double cp, double gamma) {
		int iMax = l1 + l2 + 1;
		double[] a = new double[iMax];

		int i, r, u, index;

		for (i = 0; i < iMax; i++) {
			for (r = 0; r < ((int) (Math.floor(i / 2.0) + 1.0)); r++) {
				for (u = 0; u < ((int) (Math.floor((i - 2.0 * r) / 2.0) + 1.0)); u++) {
					index = i - 2 * r - u;

					a[index] += constructATerm(i, r, u, l1, l2, pa, pb, cp,
							gamma);
				}
			}
		}

		return a;
	}

	public double[] constructGradAArray(int l1, int l2, double pa, double pb,
			double cp, double gamma) {
		int iMax = l1 + l2 + 1;
		double[] a = new double[iMax];

		int i, r, u, index;

		for (i = 0; i < iMax; i++) {
			for (r = 0; r < ((int) (Math.floor(i / 2.0) + 1.0)); r++) {
				for (u = 0; u < ((int) (Math.floor((i - 2.0 * r) / 2.0) + 1.0)); u++) {
					index = i - 2 * r - u;

					a[index] += (i - 2.0 * r + 1.0)
							/ (i - 2.0 * r - 2.0 * u + 1.0)
							* cp
							* constructATerm(i, r, u, l1, l2, pa, pb, cp, gamma);
				}
			}
		}

		return a;
	}

	/**
	 * the A term <br>
	 * <i> THO eq. 2.18 <i>
	 */
	public double constructATerm(int i, int r, int u, int l1, int l2,
			double pax, double pbx, double cpx, double gamma) {
		return (Math.pow(-1.0, i)
				* MathUtil.binomialPrefactor(i, l1, l2, pax, pbx)
				* Math.pow(-1.0, u) * MathUtil.factorial(i)
				* Math.pow(cpx, i - 2 * r - 2 * u)
				* Math.pow(0.25 / gamma, r + u) / MathUtil.factorial(r)
				/ MathUtil.factorial(u) / MathUtil.factorial(i - 2 * r - 2 * u));
	}

	/**
	 * The nuclear attraction gradient term
	 * 
	 * @param a
	 *            the coefficient of primitive Gaussian a.
	 * @param power1
	 *            the orbital powers of primitive Gaussian a.
	 * @param alpha1
	 *            the coefficient of primitive Gaussian a.
	 * @param b
	 *            the coefficient of primitive Gaussian b.
	 * @param power2
	 *            the orbital powers of primitive Gaussian b.
	 * @param alpha2
	 *            the coefficient of primitive Gaussian b.
	 * @param c
	 *            the location of nuclear centre.
	 * @return the nuclear attraction gradient.
	 */
	public Vector3D nuclearAttractionGradient(Point3D a, Power power1,
			double alpha1, Point3D b, Power power2, double alpha2, Point3D c) {
		Point3D product = IntegralsUtil.gaussianProductCenter(alpha1, a,
				alpha2, b);

		double rABSquared = a.distanceSquaredFrom(b);
		double rCPSquared = c.distanceSquaredFrom(product);

		double gamma = alpha1 + alpha2;

		double[] ax = constructAArray(power1.getL(), power2.getL(),
				product.getX() - a.getX(), product.getX() - b.getX(),
				product.getX() - c.getX(), gamma);

		double[] ay = constructAArray(power1.getM(), power2.getM(),
				product.getY() - a.getY(), product.getY() - b.getY(),
				product.getY() - c.getY(), gamma);

		double[] az = constructAArray(power1.getN(), power2.getN(),
				product.getZ() - a.getZ(), product.getZ() - b.getZ(),
				product.getZ() - c.getZ(), gamma);

		double[] dax = constructGradAArray(power1.getL(), power2.getL(),
				product.getX() - a.getX(), product.getX() - b.getX(),
				product.getX() - c.getX(), gamma);

		double[] day = constructGradAArray(power1.getM(), power2.getM(),
				product.getY() - a.getY(), product.getY() - b.getY(),
				product.getY() - c.getY(), gamma);

		double[] daz = constructGradAArray(power1.getN(), power2.getN(),
				product.getZ() - a.getZ(), product.getZ() - b.getZ(),
				product.getZ() - c.getZ(), gamma);
		double gradX = 0.0, gradY = 0.0, gradZ = 0.0;
		int i, j, k;

		for (i = 0; i < dax.length; i++) {
			for (j = 0; j < day.length; j++) {
				for (k = 0; k < daz.length; k++) {
					double fgamma = IntegralsUtil.computeFGamma(i + j + k + 1,
							rCPSquared * gamma);

					gradX += dax[i] * ay[j] * az[k] * fgamma;
					gradY += ax[i] * day[j] * az[k] * fgamma;
					gradZ += ax[i] * ay[j] * daz[k] * fgamma;
				}
			}
		}

		double factor = -4 * Math.PI
				* Math.exp(-alpha1 * alpha2 * rABSquared / gamma);

		gradX = factor * gradX;
		gradY = factor * gradY;
		gradZ = factor * gradZ;

		return new Vector3D(gradX, gradY, gradZ);
	}
}
