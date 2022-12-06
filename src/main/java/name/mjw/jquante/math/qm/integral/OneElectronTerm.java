package name.mjw.jquante.math.qm.integral;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.basis.Power;
import net.jafama.FastMath;

/**
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class OneElectronTerm implements IntegralsPackage {
	/**
	 * 1D overlap. <br>
	 * <i> Taken from http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.12 </i>
	 * 
	 * @param l1    the angular momentum number of Gaussian 1.
	 * @param l2    the angular momentum number of Gaussian 2.
	 * @param pax   the distance of Gaussian 1 to the product centre.
	 * @param pbx   the distance of Gaussian 2 to the product centre.
	 * @param gamma the sum of both Gaussian's exponent.
	 * @return the 1D overlap.
	 */
	public final double overlap1D(final int l1, final int l2, final double pax, final double pbx, final double gamma) {
		double sum = 0.0;
		final int k = 1 + (int) FastMath.floor(0.5 * (l1 + l2));

		for (int i = 0; i < k; i++) {
			sum += (MathUtil.binomialPrefactor(2 * i, l1, l2, pax, pbx) * MathUtil.factorial2(2 * i - 1))
					/ FastMath.pow(2D * gamma, i);
		}
		return sum;
	}

	/**
	 * Overlap matrix element taken from <br>
	 * <i> Taken from http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.12 </i>
	 * 
	 * @param alpha1 the coefficient of primitive Gaussian a.
	 * @param power1 the orbital powers of primitive Gaussian a.
	 * @param a      the location of primitive Gaussian a.
	 * 
	 * @param alpha2 the coefficient of primitive Gaussian b.
	 * @param power2 the orbital powers of primitive Gaussian b.
	 * @param b      the location of primitive Gaussian b.
	 * @return the Overlap integral
	 */
	public final double overlap(final double alpha1, final Power power1, final Vector3D a, final double alpha2,
			final Power power2, final Vector3D b) {

		final double radiusABSquared = a.distanceSq(b);
		final double gamma = alpha1 + alpha2;
		final Vector3D product = IntegralsUtil.gaussianProductCenter(alpha1, a, alpha2, b);

		final double wx = overlap1D(power1.l(), power2.l(), product.getX() - a.getX(), product.getX() - b.getX(),
				gamma);

		final double wy = overlap1D(power1.m(), power2.m(), product.getY() - a.getY(), product.getY() - b.getY(),
				gamma);

		final double wz = overlap1D(power1.n(), power2.n(), product.getZ() - a.getZ(), product.getZ() - b.getZ(),
				gamma);

		return (FastMath.pow(FastMath.PI / gamma, 1.5) * FastMath.exp((-alpha1 * alpha2 * radiusABSquared) / gamma) * wx
				* wy * wz);
	}

	/**
	 * The Kinetic Energy (KE) component. <br>
	 * 
	 * <i> Taken from http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.13 </i>
	 *
	 * @param alpha1 the coefficient of primitive Gaussian a.
	 * @param power1 the orbital powers of primitive Gaussian a.
	 * @param a      the location of primitive Gaussian a.
	 *
	 * @param alpha2 the coefficient of primitive Gaussian b.
	 * @param power2 the orbital powers of primitive Gaussian b.
	 * @param b      the location of primitive Gaussian b.
	 * @return the Kinetic Energy integral
	 */
	public final double kinetic(final double alpha1, final Power power1, final Vector3D a, final double alpha2,
			final Power power2, final Vector3D b) {
		final int l2 = power2.l();
		final int m2 = power2.m();
		final int n2 = power2.n();

		double term = alpha2 * (2 * (l2 + m2 + n2) + 3) * overlap(alpha1, power1, a, alpha2, power2, b);

		term += -2.0 * FastMath.pow2(alpha2)
				* (overlap(alpha1, power1, a, alpha2, new Power(l2 + 2, m2, n2), b)
						+ overlap(alpha1, power1, a, alpha2, new Power(l2, m2 + 2, n2), b)
						+ overlap(alpha1, power1, a, alpha2, new Power(l2, m2, n2 + 2), b));

		term += -0.5 * ((l2 * (l2 - 1)) * overlap(alpha1, power1, a, alpha2, new Power(l2 - 2, m2, n2), b)
				+ (m2 * (m2 - 1)) * overlap(alpha1, power1, a, alpha2, new Power(l2, m2 - 2, n2), b)
				+ (n2 * (n2 - 1)) * overlap(alpha1, power1, a, alpha2, new Power(l2, m2, n2 - 2), b));

		return term;
	}
}
