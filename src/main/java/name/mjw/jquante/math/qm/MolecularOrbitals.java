package name.mjw.jquante.math.qm;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.la.Diagonalizer;
import name.mjw.jquante.math.la.DiagonalizerFactory;

/**
 * Represents the Molecular orbitals as a coefficient matrix and the
 * corresponding eigenvalues representing the orbital energies.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class MolecularOrbitals extends Matrix {

	private static final Logger LOG = LogManager.getLogger(MolecularOrbitals.class);

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public MolecularOrbitals(int n) {
		super(n, n);
	}

	/**
	 * Get the coefficient Matrix for this MolecularOrbitals
	 * 
	 * @return the coefficient matrix
	 */
	public double[][] getCoefficients() {
		return getData();
	}

	protected double[] orbitalEnergies;

	/**
	 * Get the value of orbitalEnergies
	 * 
	 * @return the value of orbitalEnergies
	 */
	public double[] getOrbitalEnergies() {
		return orbitalEnergies;
	}

	/**
	 * Compute the MO coefficients and the orbital energies
	 * 
	 * @param hCore
	 *            the HCore matrix
	 * @param overlap
	 *            the Overlap matrix
	 */
	public void compute(HCore hCore, Overlap overlap) {
		compute((Matrix) hCore, overlap);
	}

	/**
	 * Compute the MO coefficients and the orbital energies
	 * 
	 * @param fock
	 *            the Fock matrix
	 * @param overlap
	 *            the Overlap matrix
	 */
	public void compute(Fock fock, Overlap overlap) {
		compute((Matrix) fock, overlap);
	}

	/** The actual computation is irrelevant of the type of matrix */
	private void compute(Matrix theMat, Overlap overlap) {
		Matrix x = overlap.getSHalf();
		Matrix a = x.multiply(theMat).multiply(x.transpose());

		LOG.debug("Matrix x \n" + x);
		LOG.debug("Matrix a \n" + a);

		Diagonalizer diag = DiagonalizerFactory.getInstance().getDefaultDiagonalizer();
		diag.diagonalize(a);

		LOG.debug("diag.getEigenVectors() \n" + diag.getEigenVectors());

		orbitalEnergies = diag.getEigenValues();
		this.setMatrix(diag.getEigenVectors().multiply(x).getData());
	}
}
