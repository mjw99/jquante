package name.mjw.jquante.math.qm;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.EigenDecomposition;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.SortedEigenDecomposition;

/**
 * Represents the Molecular orbitals as a coefficient matrix and the
 * corresponding eigenvalues representing the orbital energies.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class MolecularOrbitals extends Array2DRowRealMatrix {

	private static final long serialVersionUID = -3159550917810939744L;

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
		compute((RealMatrix) hCore, overlap);
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
		compute((RealMatrix) fock, overlap);
	}

	/** The actual computation is irrelevant of the type of matrix */
	private void compute(RealMatrix theMat, Overlap overlap) {
		LOG.debug("");
		RealMatrix x = overlap.getSHalf();
		LOG.debug("x: {}", x);
		RealMatrix a = x.multiply(theMat).multiply(x.transpose());
		LOG.debug("a: {}", a);

		EigenDecomposition eig = new EigenDecomposition(a);

		SortedEigenDecomposition sortedEig = new SortedEigenDecomposition(eig);

		orbitalEnergies = sortedEig.getRealEigenvalues();

		this.setSubMatrix(sortedEig.getVT().multiply(x).getData(), 0, 0);

		LOG.debug("MO values :{}", this);
	}

	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}

}
