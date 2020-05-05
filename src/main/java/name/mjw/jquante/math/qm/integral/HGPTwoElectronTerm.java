package name.mjw.jquante.math.qm.integral;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import net.jafama.FastMath;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;

import com.google.common.primitives.Doubles;

/**
 * Head-Gordon/Pople scheme of evaluating two-electron integrals.
 * 
 * The code is based upon PyQuante (<a href="http://pyquante.sf.net">
 * http://pyquante.sf.net </a>). See 'A method for two‐electron Gaussian
 * integral and integral derivative evaluation using recurrence relations'
 * <a href="http://dx.doi.org/10.1063/1.455553"> M. Head-Gordon and J. A. Pople,
 * J. Chem. Phys. <b>89</b>, 5777 (1988)</a> for more details.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class HGPTwoElectronTerm implements TwoElectronTerm {

	private final double sqrt2PI = FastMath.sqrt(2.0) * FastMath.pow(Math.PI, 1.25);

	/**
	 * 2E coulomb interactions between four contracted Gaussians.
	 */
	@Override
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c,
			ContractedGaussian d) {

		final double[] aExps = Doubles.toArray(a.getExponents());
		final double[] aCoefs = Doubles.toArray(a.getCoefficients());
		final double[] aNorms = Doubles.toArray(a.getPrimNorms());
		final int[] aPowers = { a.getPowers().getL(), a.getPowers().getM(), a.getPowers().getN() };
		final double[] aCoord = { a.getOrigin().getX(), a.getOrigin().getY(), a.getOrigin().getZ() };

		final double[] bExps = Doubles.toArray(b.getExponents());
		final double[] bCoefs = Doubles.toArray(b.getCoefficients());
		final double[] bNorms = Doubles.toArray(b.getPrimNorms());
		final int[] bPowers = { b.getPowers().getL(), b.getPowers().getM(), b.getPowers().getN() };
		final double[] bCoord = { b.getOrigin().getX(), b.getOrigin().getY(), b.getOrigin().getZ() };

		final double[] cExps = Doubles.toArray(c.getExponents());
		final double[] cCoefs = Doubles.toArray(c.getCoefficients());
		final double[] cNorms = Doubles.toArray(c.getPrimNorms());
		final int[] cPowers = { c.getPowers().getL(), c.getPowers().getM(), c.getPowers().getN() };
		final double[] cCoord = { c.getOrigin().getX(), c.getOrigin().getY(), c.getOrigin().getZ() };

		final double[] dExps = Doubles.toArray(d.getExponents());
		final double[] dCoefs = Doubles.toArray(d.getCoefficients());
		final double[] dNorms = Doubles.toArray(d.getPrimNorms());
		final int[] dPowers = { d.getPowers().getL(), d.getPowers().getM(), d.getPowers().getN() };
		final double[] dCoord = { d.getOrigin().getX(), d.getOrigin().getY(), d.getOrigin().getZ() };

		return (a.getNormalization() * b.getNormalization() * c.getNormalization() * d.getNormalization()
				* contractedHrr(aCoord, aPowers, aCoefs, aExps, aNorms, bCoord, bPowers, bCoefs, bExps, bNorms, cCoord,
						cPowers, cCoefs, cExps, cNorms, dCoord, dPowers, dCoefs, dExps, dNorms));
	}

	/**
	 * Coulomb repulsion term
	 */
	@Override
	public final double coulombRepulsion(Vector3D a, double aNorm, Power aPower, double aAlpha, Vector3D b,
			double bNorm, Power bPower, double bAlpha, Vector3D c, double cNorm, Power cPower, double cAlpha,
			Vector3D d, double dNorm, Power dPower, double dAlpha) {
		throw new UnsupportedOperationException("Not supported anymore.");
	}

	/**
	 * Contracted HRR (Horizontal Recurrence Relation)
	 * 
	 * @param a      Center of contracted Gaussian function a.
	 * @param aPower Angular momentum of Gaussian function a.
	 * @param aCoeff Coefficients of primitives in contracted Gaussian function a.
	 * @param aExps  Orbital exponents of primitives in contracted Gaussian function
	 *               a.
	 * @param aNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function a.
	 * 
	 * @param b      Center of contracted Gaussian function b.
	 * @param bPower Angular momentum of Gaussian function b.
	 * @param bCoeff Coefficients of primitives in contracted Gaussian function b.
	 * @param bExps  Orbital exponents of primitives in contracted Gaussian function
	 *               b.
	 * @param bNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function b.
	 * @param c      Center of contracted Gaussian function c.
	 * @param cPower Angular momentum of Gaussian function c.
	 * @param cCoeff Coefficients of primitives in contracted Gaussian function c.
	 * @param cExps  Orbital exponents of primitives in contracted Gaussian function
	 *               c.
	 * @param cNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function c.
	 * @param d      Center of contracted Gaussian function d.
	 * @param dPower Angular momentum of Gaussian function d.
	 * @param dCoeff Coefficients of primitives in contracted Gaussian function d.
	 * @param dExps  Orbital exponents of primitives in contracted Gaussian function
	 *               d.
	 * @param dNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function d.
	 * 
	 * @return Contribution to Horizontal Recurrence Relation.
	 */
	protected final double contractedHrr(double[] aCoord, int[] aPower, double[] aCoefs, double[] aExps,
			double[] aNorms, double[] bCoord, int[] bPower, double[] bCoefs, double[] bExps, double[] bNorms,
			double[] cCoord, int[] cPower, double[] cCoefs, double[] cExps, double[] cNorms, double[] dCoord,
			int[] dPower, double[] dCoefs, double[] dExps, double[] dNorms) {

		final int la = aPower[0];
		final int ma = aPower[1];
		final int na = aPower[2];

		final int lb = bPower[0];
		final int mb = bPower[1];
		final int nb = bPower[2];

		final int lc = cPower[0];
		final int mc = cPower[1];
		final int nc = cPower[2];

		final int ld = dPower[0];
		final int md = dPower[1];
		final int nd = dPower[2];

		if (lb > 0) {

			int[] newBPowers = { lb - 1, mb, nb };
			int[] newAPowers = { la + 1, ma, na };

			return (contractedHrr(aCoord, newAPowers, aCoefs, aExps, aNorms, bCoord, newBPowers, bCoefs, bExps, bNorms,
					cCoord, cPower, cCoefs, cExps, cNorms, dCoord, dPower, dCoefs, dExps, dNorms)

					+ (aCoord[0] - bCoord[0]) *

							contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, newBPowers, bCoefs, bExps,
									bNorms, cCoord, cPower, cCoefs, cExps, cNorms, dCoord, dPower, dCoefs, dExps,
									dNorms));

		} else if (mb > 0) {
			int[] newBPowers = { lb, mb - 1, nb };
			int[] newAPowers = { la, ma + 1, na };

			return (contractedHrr(aCoord, newAPowers, aCoefs, aExps, aNorms, bCoord, newBPowers, bCoefs, bExps, bNorms,
					cCoord, cPower, cCoefs, cExps, cNorms, dCoord, dPower, dCoefs, dExps, dNorms)

					+ (aCoord[1] - bCoord[1]) *

							contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, newBPowers, bCoefs, bExps,
									bNorms, cCoord, cPower, cCoefs, cExps, cNorms, dCoord, dPower, dCoefs, dExps,
									dNorms));

		} else if (nb > 0) {
			int[] newBPowers = { lb, mb, nb - 1 };
			int[] newAPowers = { la, ma, na + 1 };

			return (contractedHrr(aCoord, newAPowers, aCoefs, aExps, aNorms, bCoord, newBPowers, bCoefs, bExps, bNorms,
					cCoord, cPower, cCoefs, cExps, cNorms, dCoord, dPower, dCoefs, dExps, dNorms)

					+ (aCoord[2] - bCoord[2]) *

							contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, newBPowers, bCoefs, bExps,
									bNorms, cCoord, cPower, cCoefs, cExps, cNorms, dCoord, dPower, dCoefs, dExps,
									dNorms));

		} else if (ld > 0) {
			int[] newDPowers = { ld - 1, md, nd };
			int[] newCPowers = { lc + 1, mc, nc };

			return (contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bPower, bCoefs, bExps, bNorms, cCoord,
					newCPowers, cCoefs, cExps, cNorms, dCoord, newDPowers, dCoefs, dExps, dNorms)

					+ (cCoord[0] - dCoord[0]) *

							contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bPower, bCoefs, bExps, bNorms,
									cCoord, cPower, cCoefs, cExps, cNorms, dCoord, newDPowers, dCoefs, dExps, dNorms));

		} else if (md > 0) {
			int[] newDPowers = { ld, md - 1, nd };
			int[] newCPowers = { lc, mc + 1, nc };

			return (contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bPower, bCoefs, bExps, bNorms, cCoord,
					newCPowers, cCoefs, cExps, cNorms, dCoord, newDPowers, dCoefs, dExps, dNorms)

					+ (cCoord[1] - dCoord[1]) *

							contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bPower, bCoefs, bExps, bNorms,
									cCoord, cPower, cCoefs, cExps, cNorms, dCoord, newDPowers, dCoefs, dExps, dNorms));
		} else if (nd > 0) {
			int[] newDPowers = { ld, md, nd - 1 };
			int[] newCPowers = { lc, mc, nc + 1 };

			return (contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bPower, bCoefs, bExps, bNorms, cCoord,
					newCPowers, cCoefs, cExps, cNorms, dCoord, newDPowers, dCoefs, dExps, dNorms)

					+ (cCoord[2] - dCoord[2]) *

							contractedHrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bPower, bCoefs, bExps, bNorms,
									cCoord, cPower, cCoefs, cExps, cNorms, dCoord, newDPowers, dCoefs, dExps, dNorms));
		}

		return contractedVrr(aCoord, aPower, aCoefs, aExps, aNorms, bCoord, bCoefs, bExps, bNorms, cCoord, cPower,
				cCoefs, cExps, cNorms, dCoord, dCoefs, dExps, dNorms);

	}

	/**
	 * Contracted VRR (Vertical Recurrence Relation) contribution
	 * 
	 * @param a      Center of contracted Gaussian function a.
	 * @param aPower Angular momentum of Gaussian function a.
	 * @param aCoeff Coefficients of primitives in contracted Gaussian function a.
	 * @param aExps  Orbital exponents of primitives in contracted Gaussian function
	 *               a.
	 * @param aNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function a.
	 * 
	 * @param b      Center of contracted Gaussian function b. *
	 * @param bCoeff Coefficients of primitives in contracted Gaussian function b.
	 * @param bExps  Orbital exponents of primitives in contracted Gaussian function
	 *               b.
	 * @param bNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function b.
	 * 
	 * @param c      Center of contracted Gaussian function c.
	 * @param cPower Angular momentum of Gaussian function c.
	 * @param cCoeff Coefficients of primitives in contracted Gaussian function c.
	 * @param cExps  Orbital exponents of primitives in contracted Gaussian function
	 *               c.
	 * @param cNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function c.
	 * 
	 * @param d      Center of contracted Gaussian function d. *
	 * @param dCoeff Coefficients of primitives in contracted Gaussian function d.
	 * @param dExps  Orbital exponents of primitives in contracted Gaussian function
	 *               d.
	 * @param dNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function d.
	 * 
	 * @return Contribution to Vertical Recurrence Relation.
	 */
	protected final double contractedVrr(double[] aCoord, int[] aPower, double[] aCoefs, double[] aExps,
			double[] aNorms, double[] bCoord, double[] bCoefs, double[] bExps, double[] bNorms, double[] cCoord,
			int[] cPower, double[] cCoefs, double[] cExps, double[] cNorms, double[] dCoord, double[] dCoefs,
			double[] dExps, double[] dNorms) {

		double value = 0.0;

		for (int i = 0; i < aExps.length; i++) {
			for (int j = 0; j < bExps.length; j++) {
				for (int k = 0; k < cExps.length; k++) {
					for (int l = 0; l < dExps.length; l++) {

						value += aCoefs[i] * bCoefs[j] * cCoefs[k] * dCoefs[l]
								* vrr(aCoord, aNorms[i], aPower, aExps[i], bCoord, bNorms[j], bExps[j], cCoord,
										cNorms[k], cPower, cExps[k], dCoord, dNorms[l], dExps[l], 0);
					}
				}
			}
		}

		return value;
	}

	/**
	 * VRR (Vertical Recurrence Relation) contribution
	 * 
	 * @param a      Center of primitive Gaussian function a.
	 * @param aNorm  Normalisation coefficients of primitive Gaussian function a.
	 * @param aPower Angular momentum of primitive Gaussian function a.
	 * @param aAlpha Orbital exponent of primitiveGaussian function a.
	 * 
	 * @param b      Center of primitive Gaussian function b.
	 * @param bNorm  Normalisation coefficients of primitive Gaussian function b.
	 * @param bAlpha Orbital exponent of primitiveGaussian function b.
	 * 
	 * @param c      Center of primitive Gaussian function c.
	 * @param cNorm  Normalisation coefficients of primitive Gaussian function c.
	 * @param cPower Angular momentum of primitive Gaussian function c.
	 * @param cAlpha Orbital exponent of primitiveGaussian function c.
	 * 
	 * @param d      Center of primitive Gaussian function d.
	 * @param dNorm  Normalisation coefficients of primitive Gaussian function d.
	 * @param dAlpha Orbital exponent of primitiveGaussian function d.
	 * 
	 * @param m      If equal to 0, ERI is true, else greater than 0, ERI is an
	 *               auxiliary integral.
	 * @return Contribution to Vertical Recurrence Relation.
	 */
	final double vrr(double[] a, double aNorm, int[] aPower, double aAlpha, double[] b, double bNorm, double bAlpha,
			double[] c, double cNorm, int[] cPower, double cAlpha, double[] d, double dNorm, double dAlpha, int m) {

		double val;

		final double[] p = IntegralsUtil.gaussianProductCenter(aAlpha, a, bAlpha, b);
		final double[] q = IntegralsUtil.gaussianProductCenter(cAlpha, c, dAlpha, d);

		final double zeta = aAlpha + bAlpha;
		final double eta = cAlpha + dAlpha;
		final double zetaPlusEta = zeta + eta;
		final double zetaByZetaPlusEta = zeta / zetaPlusEta;
		final double etaByZetaPlusEta = eta / zetaPlusEta;

		final double[] w = IntegralsUtil.gaussianProductCenter(zeta, p, eta, q);

		final int la = aPower[0];
		final int ma = aPower[1];
		final int na = aPower[2];

		final int lc = cPower[0];
		final int mc = cPower[1];
		final int nc = cPower[2];

		if (nc > 0) {

			int[] newCPower = { lc, mc, nc - 1 };

			val = (q[2] - c[2])
					* vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m)

					+ (w[2] - q[2]) * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d,
							dNorm, dAlpha, m + 1);

			if (nc > 1) {
				int[] newCPower1 = { lc, mc, nc - 2 };

				val += 0.5 * (nc - 1) / eta
						* (vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower1, cAlpha, d, dNorm,
								dAlpha, m)
								- zetaByZetaPlusEta * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm,
										newCPower1, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			if (na > 0) {
				int[] newAPower = { la, ma, na - 1 };

				val += 0.5 * na / zetaPlusEta * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower,
						cAlpha, d, dNorm, dAlpha, m + 1);
			}
			return val;

		} else if (mc > 0) {

			int[] newCPower = { lc, mc - 1, nc };

			val = (q[1] - c[1])
					* vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m)

					+ (w[1] - q[1]) * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d,
							dNorm, dAlpha, m + 1);

			if (mc > 1) {

				int[] newCPower1 = { lc, mc - 2, nc };

				val += 0.5 * (mc - 1) / eta

						* (vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower1, cAlpha, d, dNorm,
								dAlpha, m)

								- zetaByZetaPlusEta * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm,
										newCPower1, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			if (ma > 0) {

				int[] newAPower = { la, ma - 1, na };

				val += 0.5 * ma / zetaPlusEta * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower,
						cAlpha, d, dNorm, dAlpha, m + 1);

			}
			return val;

		} else if (lc > 0) {

			int[] newCPower = { lc - 1, mc, nc };

			val = (q[0] - c[0])
					* vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d, dNorm, dAlpha, m)

					+ (w[0] - q[0]) * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower, cAlpha, d,
							dNorm, dAlpha, m + 1);

			if (lc > 1) {
				int[] newCPower1 = { lc - 2, mc, nc };

				val += 0.5 * (lc - 1) / eta
						* (vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower1, cAlpha, d, dNorm,
								dAlpha, m)

								- zetaByZetaPlusEta * vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm,
										newCPower1, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			if (la > 0) {

				int[] newAPower = { la - 1, ma, na };

				val += 0.5 * la / zetaPlusEta * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, newCPower,
						cAlpha, d, dNorm, dAlpha, m + 1);
			}
			return val;

		} else if (na > 0) {

			int[] newAPower = { la, ma, na - 1 };

			val = (p[2] - a[2])
					* vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m)

					+ (w[2] - p[2]) * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d,
							dNorm, dAlpha, m + 1);

			if (na > 1) {

				int[] newAPower1 = { la, ma, na - 2 };

				val += 0.5 * (na - 1) / zeta
						* (vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm,
								dAlpha, m)

								- etaByZetaPlusEta * vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm,
										cPower, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			return val;

		} else if (ma > 0) {
			int[] newAPower = { la, ma - 1, na };

			val = (p[1] - a[1])
					* vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m)

					+ (w[1] - p[1]) * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d,
							dNorm, dAlpha, m + 1);

			if (ma > 1) {
				int[] newAPower1 = { la, ma - 2, na };

				val += 0.5 * (ma - 1) / zeta
						* (vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm,
								dAlpha, m)

								- etaByZetaPlusEta * vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm,
										cPower, cAlpha, d, dNorm, dAlpha, m + 1));
			}

			return val;

		} else if (la > 0) {

			int[] newAPower = { la - 1, ma, na };

			val = (p[0] - a[0])
					* vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m)

					+ (w[0] - p[0]) * vrr(a, aNorm, newAPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d,
							dNorm, dAlpha, m + 1);

			if (la > 1) {
				int[] newAPower1 = { la - 2, ma, na };

				val += 0.5 * (la - 1) / zeta
						* (vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm,
								dAlpha, m)

								- etaByZetaPlusEta * vrr(a, aNorm, newAPower1, aAlpha, b, bNorm, bAlpha, c, cNorm,
										cPower, cAlpha, d, dNorm, dAlpha, m + 1));

			}

			return val;
		}

		final double rab2 = MathUtil.distanceSq(a, b);
		final double Kab = sqrt2PI / zeta * FastMath.exp(-aAlpha * bAlpha / zeta * rab2);

		final double rcd2 = MathUtil.distanceSq(c, d);
		final double Kcd = sqrt2PI / eta * FastMath.exp(-cAlpha * dAlpha / eta * rcd2);

		final double rpq2 = MathUtil.distanceSq(p, q);
		final double T = zeta * eta / zetaPlusEta * rpq2;

		return aNorm * bNorm * cNorm * dNorm * Kab * Kcd / FastMath.sqrt(zetaPlusEta)
				* IntegralsUtil.computeFGamma(m, T);
	}

	@Override
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c, ContractedGaussian d,
			Density density, RealMatrix jMat, RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public final double coulombRepulsion(double[] aCoord, double aNorm, int[] aPowers, double aAlpha, double[] bCoord,
			double bNorm, int[] bPowers, double bAlpha, double[] cCoord, double cNorm, int[] cPowers, double cAlpha,
			double[] dCoord, double dNorm, int[] dPowers, double dAlpha) {

		return vrr(aCoord, aNorm, aPowers, aAlpha, bCoord, bNorm, bAlpha, cCoord, cNorm, cPowers, cAlpha, dCoord, dNorm,
				dAlpha, 0);

	}
}
