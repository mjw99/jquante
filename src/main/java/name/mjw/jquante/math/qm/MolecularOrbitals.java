package name.mjw.jquante.math.qm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Represents the Molecular orbitals as a coefficient matrix and the
 * corresponding eigenvalues representing the orbital energies.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class MolecularOrbitals extends Array2DRowRealMatrix {

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
		LOG.debug("x: " + x);
		RealMatrix a = x.multiply(theMat).multiply(x.transpose());
		LOG.debug("a: " + a);

		// Ensure the matrix is symmetric.
		// This is so that the EigenDecomposition.findEigenVectors() execution is followed
		// and the resulting eigenvectors are sorted by size.
		final double tol = 10 * a.getRowDimension() * a.getColumnDimension() * Precision.EPSILON;
		for (int i = 0; i < a.getRowDimension(); i++) {
			for (int j = 0; j < a.getColumnDimension(); j++) {
				final double mij = a.getEntry(i, j);
				final double mji = a.getEntry(j, i);
				if (FastMath.abs(mij - mji) > FastMath.max(FastMath.abs(mij), FastMath.abs(mji)) * tol) {
					a.setEntry(i, j, a.getEntry(j, i));
				}
			}
		}

		LOG.debug("is a symmetric? :" + MatrixUtils.isSymmetric(a, tol));

		EigenDecomposition eig = new EigenDecomposition(a);

		LOG.debug("eig.getD(): " + eig.getD());
		LOG.debug("eig.getVT(): " + eig.getVT());

		orbitalEnergies = eig.getRealEigenvalues();

		int j = eig.getRealEigenvalues().length - 1;
		RealMatrix sortedVT = MatrixUtils.createRealMatrix(orbitalEnergies.length, orbitalEnergies.length);

		// Reverse the order of both eigenvalues and eigenvectors
		for (int i = 0; i < eig.getRealEigenvalues().length; i++) {
			orbitalEnergies[i] = eig.getRealEigenvalue(j);
			sortedVT.setRowVector(i, eig.getVT().getRowVector(j));
			j--;
		}

		this.setSubMatrix(sortedVT.multiply(x).getData(), 0, 0);

		LOG.debug("MO::values :" + this);
	}

}
