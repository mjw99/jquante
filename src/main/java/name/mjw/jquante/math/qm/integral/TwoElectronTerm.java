/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package name.mjw.jquante.math.qm.integral;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.math.qm.Density;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

/**
 * Top level interface for evaluating a 2E-integral term.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public abstract class TwoElectronTerm implements IntegralsPackage {

	/**
	 * 2E coulomb interactions between 4 contracted Gaussians
	 */
	public abstract double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d);

	/** The column repulsion term between four centered Gaussians */
	public abstract double coulombRepulsion(Point3D a, double aNorm,
			Power aPower, double aAlpha, Point3D b, double bNorm, Power bPower,
			double bAlpha, Point3D c, double cNorm, Power cPower,
			double cAlpha, Point3D d, double dNorm, Power dPower, double dAlpha);

	/**
	 * 2E coulomb interactions between 4 contracted Gaussians
	 */
	public abstract double coulomb(ContractedGaussian a, ContractedGaussian b,
			ContractedGaussian c, ContractedGaussian d, Density density,
			Matrix jMat, Matrix kMat);

}
