package name.mjw.jquante.math.qm.integral;

import java.util.ArrayList;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.RealMatrix;

import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

/**
 * Two electron integrals using RYS quadrature. The code is based upon PyQuante
 * (<a href="http://pyquante.sf.net"> http://pyquante.sf.net </a>).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class RysTwoElectronTerm extends TwoElectronTerm {
	/**
	 * 2E coulomb interactions between 4 contracted Gaussians
	 */
	@Override
	public double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d) {
		ArrayList<Double> aExps, aCoefs, aNorms, bExps, bCoefs, bNorms, cExps, cCoefs, cNorms, dExps, dCoefs, dNorms;

		aExps = a.getExponents();
		aCoefs = a.getCoefficients();
		aNorms = a.getPrimNorms();
		Vector3D aOrigin = a.getOrigin();
		Power aPower = a.getPowers();

		bExps = b.getExponents();
		bCoefs = b.getCoefficients();
		bNorms = b.getPrimNorms();
		Vector3D bOrigin = b.getOrigin();
		Power bPower = b.getPowers();

		cExps = c.getExponents();
		cCoefs = c.getCoefficients();
		cNorms = c.getPrimNorms();
		Vector3D cOrigin = c.getOrigin();
		Power cPower = c.getPowers();

		dExps = d.getExponents();
		dCoefs = d.getCoefficients();
		dNorms = d.getPrimNorms();
		Vector3D dOrigin = d.getOrigin();
		Power dPower = d.getPowers();

		// TODO:
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * coulomb repulsion term
	 */
	@Override
	public double coulombRepulsion(Vector3D a, double aNorm, Power aPower,
			double aAlpha, Vector3D b, double bNorm, Power bPower,
			double bAlpha, Vector3D c, double cNorm, Power cPower,
			double cAlpha, Vector3D d, double dNorm, Power dPower, double dAlpha) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d, Density density,
			RealMatrix jMat, RealMatrix kMat) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
