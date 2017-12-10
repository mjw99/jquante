package name.mjw.jquante.math.qm.integral;

import java.util.ArrayList;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

/**
 * The Huzinaga integral package.
 * 
 * The equations herein are based upon: <br>
 * 'Gaussian Expansion Methods for Molecular Orbitals.' <a
 * href="http://dx.doi.org/10.1143/JPSJ.21.2313"> H. Taketa, S. Huzinaga, and K.
 * O-ohata. <i> H. Phys. Soc. Japan, </i> <b>21</b>, 2313, 1966.)</a> [THO
 * paper]. <br>
 * and PyQuante (<a href="http://pyquante.sf.net"> http://pyquante.sf.net </a>).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class HuzinagaTwoElectronTerm extends TwoElectronTerm {
	/**
	 * 2E coulomb interactions between four contracted Gaussians
	 */
	@Override
	public double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d) {

		double jij = 0.0;

		int i;
		int j;
		int k;
		int l;
		double iaExp;
		double iaCoef;
		double iaNorm;
		double jbExp;
		double jbCoef;
		double jbNorm;
		double kcExp;
		double kcCoef;
		double kcNorm;
		double repulsionTerm;

		ArrayList<Double> aExps = a.getExponents();
		ArrayList<Double> aCoefs = a.getCoefficients();
		ArrayList<Double> aNorms = a.getPrimNorms();
		Point3D aOrigin = a.getOrigin();
		Power aPower = a.getPowers();

		ArrayList<Double> bExps = b.getExponents();
		ArrayList<Double> bCoefs = b.getCoefficients();
		ArrayList<Double> bNorms = b.getPrimNorms();
		Point3D bOrigin = b.getOrigin();
		Power bPower = b.getPowers();

		ArrayList<Double> cExps = c.getExponents();
		ArrayList<Double> cCoefs = c.getCoefficients();
		ArrayList<Double> cNorms = c.getPrimNorms();
		Point3D cOrigin = c.getOrigin();
		Power cPower = c.getPowers();

		ArrayList<Double> dExps = d.getExponents();
		ArrayList<Double> dCoefs = d.getCoefficients();
		ArrayList<Double> dNorms = d.getPrimNorms();
		Point3D dOrigin = d.getOrigin();
		Power dPower = d.getPowers();

		int asz = aExps.size();
		int bsz = bExps.size();
		int csz = cExps.size();
		int dsz = dExps.size();

		for (i = 0; i < asz; i++) {
			iaCoef = aCoefs.get(i);
			iaExp = aExps.get(i);
			iaNorm = aNorms.get(i);

			for (j = 0; j < bsz; j++) {
				jbCoef = bCoefs.get(j);
				jbExp = bExps.get(j);
				jbNorm = bNorms.get(j);

				for (k = 0; k < csz; k++) {
					kcCoef = cCoefs.get(k);
					kcExp = cExps.get(k);
					kcNorm = cNorms.get(k);

					for (l = 0; l < dsz; l++) {
						repulsionTerm = coulombRepulsion(aOrigin, iaNorm,
								aPower, iaExp, bOrigin, jbNorm, bPower, jbExp,
								cOrigin, kcNorm, cPower, kcExp, dOrigin,
								dNorms.get(l), dPower, dExps.get(l));

						jij += iaCoef * jbCoef * kcCoef * dCoefs.get(l)
								* repulsionTerm;
					}
				}
			}
		}

		return (a.getNormalization() * b.getNormalization()
				* c.getNormalization() * d.getNormalization() * jij);
	}

	/**
	 * coulomb repulsion term
	 */
	@Override
	public double coulombRepulsion(Point3D a, double aNorm, Power aPower,
			double aAlpha, Point3D b, double bNorm, Power bPower,
			double bAlpha, Point3D c, double cNorm, Power cPower,
			double cAlpha, Point3D d, double dNorm, Power dPower, double dAlpha) {

		double sum = 0.0;
		int i;
		int j;
		int k;

		double radiusABSquared = a.distanceSquaredFrom(b);
		double radiusCDSquared = c.distanceSquaredFrom(d);

		Point3D p = IntegralsUtil.gaussianProductCenter(aAlpha, a, bAlpha, b);
		Point3D q = IntegralsUtil.gaussianProductCenter(cAlpha, c, dAlpha, d);

		double radiusPQSquared = p.distanceSquaredFrom(q);

		double gamma1 = aAlpha + bAlpha;
		double gamma2 = cAlpha + dAlpha;
		double delta = 0.25 * (1 / gamma1 + 1 / gamma2);

		double[] bx = constructBArray(aPower.getL(), bPower.getL(),
				cPower.getL(), dPower.getL(), p.getX(), a.getX(), b.getX(),
				q.getX(), c.getX(), d.getX(), gamma1, gamma2, delta);

		double[] by = constructBArray(aPower.getM(), bPower.getM(),
				cPower.getM(), dPower.getM(), p.getY(), a.getY(), b.getY(),
				q.getY(), c.getY(), d.getY(), gamma1, gamma2, delta);

		double[] bz = constructBArray(aPower.getN(), bPower.getN(),
				cPower.getN(), dPower.getN(), p.getZ(), a.getZ(), b.getZ(),
				q.getZ(), c.getZ(), d.getZ(), gamma1, gamma2, delta);

		for (i = 0; i < bx.length; i++) {
			for (j = 0; j < by.length; j++) {
				for (k = 0; k < bz.length; k++) {
					sum += bx[i]
							* by[j]
							* bz[k]
							* IntegralsUtil.computeFGamma(i + j + k, 0.25
									* radiusPQSquared / delta);
				}
			}
		}

		return (2 * Math.pow(Math.PI, 2.5)
				/ (gamma1 * gamma2 * Math.sqrt(gamma1 + gamma2))
				* Math.exp(-aAlpha * bAlpha * radiusABSquared / gamma1)
				* Math.exp(-cAlpha * dAlpha * radiusCDSquared / gamma2) * sum
				* aNorm * bNorm * cNorm * dNorm);
	}

	/**
	 * Construct B array.
	 * 
	 * <i> http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.22 </i>
	 */
	private double[] constructBArray(int l1, int l2, int l3, int l4, double p,
			double a, double b, double q, double c, double d, double g1,
			double g2, double delta) {

		int i1;
		int i2;
		int r1;
		int r2;
		int u;
		int index;

		int iMax = l1 + l2 + l3 + l4 + 1;
		double[] bArr = new double[iMax];

		for (i1 = 0; i1 < (l1 + l2 + 1); i1++) {
			for (i2 = 0; i2 < (l3 + l4 + 1); i2++) {
				for (r1 = 0; r1 < (i1 / 2 + 1); r1++) {
					for (r2 = 0; r2 < (i2 / 2 + 1); r2++) {
						for (u = 0; u < ((i1 + i2) / 2 - r1 - r2 + 1); u++) {
							index = i1 + i2 - 2 * (r1 + r2) - u;

							bArr[index] += constructBTerm(i1, i2, r1, r2, u,
									l1, l2, l3, l4, p, a, b, q, c, d, g1, g2,
									delta);
						}
					}
				}
			}
		}

		return bArr;
	}

	/**
	 * Construct the B term
	 * 
	 * <i> http://dx.doi.org/10.1143/JPSJ.21.2313 eq. 2.22 </i>
	 */
	private double constructBTerm(int i1, int i2, int r1, int r2, int u,
			int l1, int l2, int l3, int l4, double px, double ax, double bx,
			double qx, double cx, double dx, double gamma1, double gamma2,
			double delta) {

		return (functionB(i1, l1, l2, px, ax, bx, r1, gamma1)
				* Math.pow(-1, i2)
				* functionB(i2, l3, l4, qx, cx, dx, r2, gamma2)
				* Math.pow(-1, u)
				* MathUtil.factorialRatioSquared(i1 + i2 - 2 * (r1 + r2), u)
				* Math.pow(qx - px, i1 + i2 - 2 * (r1 + r2) - 2 * u) / Math
					.pow(delta, i1 + i2 - 2 * (r1 + r2) - u));
	}

	/**
	 * the function B, taken from PyQuante
	 */
	private double functionB(int i, int l1, int l2, double p, double a,
			double b, int r, double g) {
		return (MathUtil.binomialPrefactor(i, l1, l2, p - a, p - b) * functionB0(
				i, r, g));
	}

	/**
	 * the function B0, taken from PyQuante
	 */
	private double functionB0(int i, int r, double g) {
		return (MathUtil.factorialRatioSquared(i, r) * Math.pow(4 * g, r - i));
	}

	@Override
	public double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d, Density density,
			Matrix jMat, Matrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
