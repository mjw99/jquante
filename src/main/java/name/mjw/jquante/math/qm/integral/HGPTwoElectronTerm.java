package name.mjw.jquante.math.qm.integral;

import java.util.ArrayList;

import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import net.jafama.FastMath;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealMatrix;

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

	/** Creates a new instance of HGPTwoElectronTerm. */
	public HGPTwoElectronTerm() {
	}

	/** Precomputed constant: sqrt(2) * PI^(5/4), used in the HGP 2E integral prefactor. */
	private final double sqrt2PI = FastMath.sqrt(2.0) * FastMath.pow(Math.PI, 1.25);

	/**
	 * 2E coulomb interactions between four contracted Gaussians using the
	 * Head-Gordon/Pople scheme.
	 *
	 * @param a the first contracted Gaussian.
	 * @param b the second contracted Gaussian.
	 * @param c the third contracted Gaussian.
	 * @param d the fourth contracted Gaussian.
	 * @return the two-electron Coulomb repulsion integral (ab|cd).
	 */
	@Override
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c,
			ContractedGaussian d) {
		return (a.getNormalization() * b.getNormalization() * c.getNormalization() * d.getNormalization()
				* contractedHrr(a.getOrigin(), a.getPowers(), a.getCoefficients(), a.getExponents(), a.getPrimNorms(),
						b.getOrigin(), b.getPowers(), b.getCoefficients(), b.getExponents(), b.getPrimNorms(),
						c.getOrigin(), c.getPowers(), c.getCoefficients(), c.getExponents(), c.getPrimNorms(),
						d.getOrigin(), d.getPowers(), d.getCoefficients(), d.getExponents(), d.getPrimNorms()));
	}

	/**
	 * Coulomb repulsion term between four primitive Gaussians using the
	 * Head-Gordon/Pople VRR scheme.
	 *
	 * @param a      the center of primitive Gaussian a.
	 * @param aNorm  the normalisation factor of primitive Gaussian a.
	 * @param aPower the angular momentum powers of primitive Gaussian a.
	 * @param aAlpha the exponent of primitive Gaussian a.
	 * @param b      the center of primitive Gaussian b.
	 * @param bNorm  the normalisation factor of primitive Gaussian b.
	 * @param bPower the angular momentum powers of primitive Gaussian b.
	 * @param bAlpha the exponent of primitive Gaussian b.
	 * @param c      the center of primitive Gaussian c.
	 * @param cNorm  the normalisation factor of primitive Gaussian c.
	 * @param cPower the angular momentum powers of primitive Gaussian c.
	 * @param cAlpha the exponent of primitive Gaussian c.
	 * @param d      the center of primitive Gaussian d.
	 * @param dNorm  the normalisation factor of primitive Gaussian d.
	 * @param dPower the angular momentum powers of primitive Gaussian d.
	 * @param dAlpha the exponent of primitive Gaussian d.
	 * @return the two-electron Coulomb repulsion integral value.
	 */
	@Override
	public final double coulombRepulsion(Vector3D a, double aNorm, Power aPower, double aAlpha, Vector3D b,
			double bNorm, Power bPower, double bAlpha, Vector3D c, double cNorm, Power cPower, double cAlpha,
			Vector3D d, double dNorm, Power dPower, double dAlpha) {
		return vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, 0);
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
	protected final double contractedHrr(Vector3D a, Power aPower, ArrayList<Double> aCoeff, ArrayList<Double> aExps,
			ArrayList<Double> aNorms, Vector3D b, Power bPower, ArrayList<Double> bCoeff, ArrayList<Double> bExps,
			ArrayList<Double> bNorms, Vector3D c, Power cPower, ArrayList<Double> cCoeff, ArrayList<Double> cExps,
			ArrayList<Double> cNorms, Vector3D d, Power dPower, ArrayList<Double> dCoeff, ArrayList<Double> dExps,
			ArrayList<Double> dNorms) {

		final int la = aPower.l();
		final int ma = aPower.m();
		final int na = aPower.n();

		final int lb = bPower.l();
		final int mb = bPower.m();
		final int nb = bPower.n();

		final int lc = cPower.l();
		final int mc = cPower.m();
		final int nc = cPower.n();

		final int ld = dPower.l();
		final int md = dPower.m();
		final int nd = dPower.n();

		if (lb > 0) {
			Power newBPower = new Power(lb - 1, mb, nb);
			return (contractedHrr(a, new Power(la + 1, ma, na), aCoeff, aExps, aNorms, b, newBPower, bCoeff, bExps,
					bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms)
					+ (a.getX() - b.getX()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, newBPower, bCoeff,
							bExps, bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms));
		} else if (mb > 0) {
			Power newBPower = new Power(lb, mb - 1, nb);
			return (contractedHrr(a, new Power(la, ma + 1, na), aCoeff, aExps, aNorms, b, newBPower, bCoeff, bExps,
					bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms)
					+ (a.getY() - b.getY()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, newBPower, bCoeff,
							bExps, bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms));
		} else if (nb > 0) {
			Power newBPower = new Power(lb, mb, nb - 1);
			return (contractedHrr(a, new Power(la, ma, na + 1), aCoeff, aExps, aNorms, b, newBPower, bCoeff, bExps,
					bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms)
					+ (a.getZ() - b.getZ()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, newBPower, bCoeff,
							bExps, bNorms, c, cPower, cCoeff, cExps, cNorms, d, dPower, dCoeff, dExps, dNorms));
		} else if (ld > 0) {
			Power newDPower = new Power(ld - 1, md, nd);
			return (contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps, bNorms, c,
					new Power(lc + 1, mc, nc), cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms)
					+ (c.getX() - d.getX()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps,
							bNorms, c, cPower, cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms));
		} else if (md > 0) {
			Power newDPower = new Power(ld, md - 1, nd);
			return (contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps, bNorms, c,
					new Power(lc, mc + 1, nc), cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms)
					+ (c.getY() - d.getY()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps,
							bNorms, c, cPower, cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms));
		} else if (nd > 0) {
			Power newDPower = new Power(ld, md, nd - 1);
			return (contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps, bNorms, c,
					new Power(lc, mc, nc + 1), cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms)
					+ (c.getZ() - d.getZ()) * contractedHrr(a, aPower, aCoeff, aExps, aNorms, b, bPower, bCoeff, bExps,
							bNorms, c, cPower, cCoeff, cExps, cNorms, d, newDPower, dCoeff, dExps, dNorms));
		}

		return contractedVrr(a, aPower, aCoeff, aExps, aNorms, b, bCoeff, bExps, bNorms, c, cPower, cCoeff, cExps,
				cNorms, d, dCoeff, dExps, dNorms);
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
	 * @param b      Center of contracted Gaussian function b.
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
	 * @param d      Center of contracted Gaussian function d.
	 * @param dCoeff Coefficients of primitives in contracted Gaussian function d.
	 * @param dExps  Orbital exponents of primitives in contracted Gaussian function
	 *               d.
	 * @param dNorms Primitive normalisation coefficients of primitives in
	 *               contracted Gaussian function d.
	 * 
	 * @return Contribution to Vertical Recurrence Relation.
	 */
	protected final double contractedVrr(Vector3D a, Power aPower, ArrayList<Double> aCoeff, ArrayList<Double> aExps,
			ArrayList<Double> aNorms, Vector3D b, ArrayList<Double> bCoeff, ArrayList<Double> bExps,
			ArrayList<Double> bNorms, Vector3D c, Power cPower, ArrayList<Double> cCoeff, ArrayList<Double> cExps,
			ArrayList<Double> cNorms, Vector3D d, ArrayList<Double> dCoeff, ArrayList<Double> dExps,
			ArrayList<Double> dNorms) {

		double value = 0.0;

		double iaExp;
		double iaCoef;
		double iaNorm;

		double jbExp;
		double jbCoef;
		double jbNorm;

		double kcExp;
		double kcCoef;
		double kcNorm;

		for (int i = 0; i < aExps.size(); i++) {
			iaCoef = aCoeff.get(i);
			iaExp = aExps.get(i);
			iaNorm = aNorms.get(i);

			for (int j = 0; j < bExps.size(); j++) {
				jbCoef = bCoeff.get(j);
				jbExp = bExps.get(j);
				jbNorm = bNorms.get(j);

				for (int k = 0; k < cExps.size(); k++) {
					kcCoef = cCoeff.get(k);
					kcExp = cExps.get(k);
					kcNorm = cNorms.get(k);

					for (int l = 0; l < dExps.size(); l++) {
						value += iaCoef * jbCoef * kcCoef * dCoeff.get(l) * vrr(a, iaNorm, aPower, iaExp, b,
								jbNorm, jbExp, c, kcNorm, cPower, kcExp, d, dNorms.get(l), dExps.get(l), 0);
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
	final double vrr(Vector3D a, double aNorm, Power aPower, double aAlpha, Vector3D b, double bNorm,
			double bAlpha, Vector3D c, double cNorm, Power cPower, double cAlpha, Vector3D d, double dNorm,
			double dAlpha, int m) {

		// All quantities below are constant for this primitive quartet.
		// Precompute them once here so that vrrRecurse() never recomputes them.
		final double zeta = aAlpha + bAlpha;
		final double eta  = cAlpha + dAlpha;
		final double zetaPlusEta = zeta + eta;

		final double ax = a.getX(), ay = a.getY(), az = a.getZ();
		final double cx = c.getX(), cy = c.getY(), cz = c.getZ();

		final double px = (aAlpha * ax + bAlpha * b.getX()) / zeta;
		final double py = (aAlpha * ay + bAlpha * b.getY()) / zeta;
		final double pz = (aAlpha * az + bAlpha * b.getZ()) / zeta;

		final double qx = (cAlpha * cx + dAlpha * d.getX()) / eta;
		final double qy = (cAlpha * cy + dAlpha * d.getY()) / eta;
		final double qz = (cAlpha * cz + dAlpha * d.getZ()) / eta;

		final double wx = (zeta * px + eta * qx) / zetaPlusEta;
		final double wy = (zeta * py + eta * qy) / zetaPlusEta;
		final double wz = (zeta * pz + eta * qz) / zetaPlusEta;

		// Geometric differences used by the VRR recurrence (ket side uses q−c and w−q;
		// bra side uses p−a and w−p).
		final double pax = px - ax, pay = py - ay, paz = pz - az;
		final double wpx = wx - px, wpy = wy - py, wpz = wz - pz;
		final double qcx = qx - cx, qcy = qy - cy, qcz = qz - cz;
		final double wqx = wx - qx, wqy = wy - qy, wqz = wz - qz;

		// Recurrence prefactors.
		final double oneover2zeta        = 0.5 / zeta;
		final double oneover2eta         = 0.5 / eta;
		final double oneover2zetaPlusEta = 0.5 / zetaPlusEta;
		final double zetaByZetaPlusEta   = zeta / zetaPlusEta;
		final double etaByZetaPlusEta    = eta  / zetaPlusEta;

		// Base-case prefactor and Boys-function argument — both constant for this quartet.
		final double rab2 = a.distanceSq(b);
		final double rcd2 = c.distanceSq(d);
		final double Kab  = sqrt2PI / zeta * FastMath.exp(-aAlpha * bAlpha / zeta * rab2);
		final double Kcd  = sqrt2PI / eta  * FastMath.exp(-cAlpha * dAlpha / eta  * rcd2);
		final double rpq2 = (px - qx) * (px - qx) + (py - qy) * (py - qy) + (pz - qz) * (pz - qz);
		final double T    = zeta * eta / zetaPlusEta * rpq2;
		final double preScaled = aNorm * bNorm * cNorm * dNorm * Kab * Kcd / FastMath.sqrt(zetaPlusEta);

		return vrrRecurse(
				aPower.l(), aPower.m(), aPower.n(),
				cPower.l(), cPower.m(), cPower.n(),
				m,
				pax, pay, paz, wpx, wpy, wpz,
				qcx, qcy, qcz, wqx, wqy, wqz,
				oneover2zeta, oneover2eta, oneover2zetaPlusEta,
				zetaByZetaPlusEta, etaByZetaPlusEta,
				preScaled, T);
	}

	/**
	 * Recursive core of the VRR (Vertical Recurrence Relation).
	 *
	 * <p>All quartet-constant geometry and prefactor quantities are passed in as
	 * plain {@code double} scalars so they are computed exactly once per primitive
	 * quartet by {@link #vrr} and never recomputed during the recursion.
	 *
	 * @param la angular momentum x of bra center a
	 * @param ma angular momentum y of bra center a
	 * @param na angular momentum z of bra center a
	 * @param lc angular momentum x of ket center c
	 * @param mc angular momentum y of ket center c
	 * @param nc angular momentum z of ket center c
	 * @param m  auxiliary integral order
	 * @param pax p.x − a.x
	 * @param pay p.y − a.y
	 * @param paz p.z − a.z
	 * @param wpx w.x − p.x
	 * @param wpy w.y − p.y
	 * @param wpz w.z − p.z
	 * @param qcx q.x − c.x
	 * @param qcy q.y − c.y
	 * @param qcz q.z − c.z
	 * @param wqx w.x − q.x
	 * @param wqy w.y − q.y
	 * @param wqz w.z − q.z
	 * @param oneover2zeta        1 / (2 * zeta)
	 * @param oneover2eta         1 / (2 * eta)
	 * @param oneover2zetaPlusEta 1 / (2 * (zeta + eta))
	 * @param zetaByZetaPlusEta   zeta / (zeta + eta)
	 * @param etaByZetaPlusEta    eta  / (zeta + eta)
	 * @param preScaled           aNorm * bNorm * cNorm * dNorm * Kab * Kcd / sqrt(zeta+eta)
	 * @param T                   Boys-function argument
	 * @return VRR integral value for the given angular momenta and order m
	 */
	private double vrrRecurse(
			int la, int ma, int na,
			int lc, int mc, int nc,
			int m,
			double pax, double pay, double paz,
			double wpx, double wpy, double wpz,
			double qcx, double qcy, double qcz,
			double wqx, double wqy, double wqz,
			double oneover2zeta,
			double oneover2eta,
			double oneover2zetaPlusEta,
			double zetaByZetaPlusEta,
			double etaByZetaPlusEta,
			double preScaled,
			double T) {
		double val;

		if (nc > 0) {
			val = qcz * vrrRecurse(la, ma, na, lc, mc, nc - 1, m,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T)
				+ wqz * vrrRecurse(la, ma, na, lc, mc, nc - 1, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			if (nc > 1) {
				val += oneover2eta * (nc - 1) * (
						vrrRecurse(la, ma, na, lc, mc, nc - 2, m,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T)
						- zetaByZetaPlusEta * vrrRecurse(la, ma, na, lc, mc, nc - 2, m + 1,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T));
			}
			if (na > 0) {
				val += oneover2zetaPlusEta * na * vrrRecurse(la, ma, na - 1, lc, mc, nc - 1, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			}
			return val;
		} else if (mc > 0) {
			val = qcy * vrrRecurse(la, ma, na, lc, mc - 1, nc, m,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T)
				+ wqy * vrrRecurse(la, ma, na, lc, mc - 1, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			if (mc > 1) {
				val += oneover2eta * (mc - 1) * (
						vrrRecurse(la, ma, na, lc, mc - 2, nc, m,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T)
						- zetaByZetaPlusEta * vrrRecurse(la, ma, na, lc, mc - 2, nc, m + 1,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T));
			}
			if (ma > 0) {
				val += oneover2zetaPlusEta * ma * vrrRecurse(la, ma - 1, na, lc, mc - 1, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			}
			return val;
		} else if (lc > 0) {
			val = qcx * vrrRecurse(la, ma, na, lc - 1, mc, nc, m,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T)
				+ wqx * vrrRecurse(la, ma, na, lc - 1, mc, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			if (lc > 1) {
				val += oneover2eta * (lc - 1) * (
						vrrRecurse(la, ma, na, lc - 2, mc, nc, m,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T)
						- zetaByZetaPlusEta * vrrRecurse(la, ma, na, lc - 2, mc, nc, m + 1,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T));
			}
			if (la > 0) {
				val += oneover2zetaPlusEta * la * vrrRecurse(la - 1, ma, na, lc - 1, mc, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			}
			return val;
		} else if (na > 0) {
			val = paz * vrrRecurse(la, ma, na - 1, lc, mc, nc, m,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T)
				+ wpz * vrrRecurse(la, ma, na - 1, lc, mc, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			if (na > 1) {
				val += oneover2zeta * (na - 1) * (
						vrrRecurse(la, ma, na - 2, lc, mc, nc, m,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T)
						- etaByZetaPlusEta * vrrRecurse(la, ma, na - 2, lc, mc, nc, m + 1,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T));
			}
			return val;
		} else if (ma > 0) {
			val = pay * vrrRecurse(la, ma - 1, na, lc, mc, nc, m,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T)
				+ wpy * vrrRecurse(la, ma - 1, na, lc, mc, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			if (ma > 1) {
				val += oneover2zeta * (ma - 1) * (
						vrrRecurse(la, ma - 2, na, lc, mc, nc, m,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T)
						- etaByZetaPlusEta * vrrRecurse(la, ma - 2, na, lc, mc, nc, m + 1,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T));
			}
			return val;
		} else if (la > 0) {
			val = pax * vrrRecurse(la - 1, ma, na, lc, mc, nc, m,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T)
				+ wpx * vrrRecurse(la - 1, ma, na, lc, mc, nc, m + 1,
						pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
						oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
						preScaled, T);
			if (la > 1) {
				val += oneover2zeta * (la - 1) * (
						vrrRecurse(la - 2, ma, na, lc, mc, nc, m,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T)
						- etaByZetaPlusEta * vrrRecurse(la - 2, ma, na, lc, mc, nc, m + 1,
								pax, pay, paz, wpx, wpy, wpz, qcx, qcy, qcz, wqx, wqy, wqz,
								oneover2zeta, oneover2eta, oneover2zetaPlusEta, zetaByZetaPlusEta, etaByZetaPlusEta,
								preScaled, T));
			}
			return val;
		}

		// Base case: all angular momenta are zero.
		return preScaled * IntegralsUtil.computeFGamma(m, T);
	}

	/**
	 * Not yet implemented. Throws {@link UnsupportedOperationException}.
	 *
	 * @param a       the first contracted Gaussian.
	 * @param b       the second contracted Gaussian.
	 * @param c       the third contracted Gaussian.
	 * @param d       the fourth contracted Gaussian.
	 * @param density the current density matrix (unused).
	 * @param jMat    the Coulomb matrix to accumulate into (unused).
	 * @param kMat    the exchange matrix to accumulate into (unused).
	 * @return never returns normally.
	 * @throws UnsupportedOperationException always.
	 */
	@Override
	public final double coulomb(ContractedGaussian a, ContractedGaussian b, ContractedGaussian c, ContractedGaussian d,
			Density density, RealMatrix jMat, RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
