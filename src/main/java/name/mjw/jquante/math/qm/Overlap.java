package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.DefaultRealMatrixChangingVisitor;
import org.hipparchus.linear.EigenDecompositionSymmetric;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.MathUtil;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import net.jafama.FastMath;

/**
 * Represents the Overlap matrix (S).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Overlap extends Array2DRowRealMatrix {

	/** Logger object. */
	private static final Logger LOG = LogManager.getLogger(Overlap.class);

	/** Eclipse-generated serialVersionUID. */
	private static final long serialVersionUID = 5209272241805533800L;

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public Overlap(int n) {
		super(n, n);
	}

	/** Cached S^(-1/2) symmetric orthogonalisation matrix. */
	private transient RealMatrix sHalf = null;

	/**
	 * Get the S^1/2 matrix
	 *
	 * Symmetric orthogonalization of the real symmetric matrix X (this). This is
	 * given by <code>U(1/sqrt(lambda))U'</code>, where lambda, U are the
	 * eigenvalues, vectors.
	 *
	 *
	 * @return return the symmetric orthogonalization matrix (S half)
	 */
	public RealMatrix getSHalf() {
		if (sHalf == null) {

			LOG.debug("Overlap::this " + this);
			// Copy and zero near-zero entries to prevent NonSymmetricMatrixException
			// when floating-point asymmetry between S[i][j] and S[j][i] exceeds
			// the default tolerance. See https://github.com/Hipparchus-Math/hipparchus/issues/365
			RealMatrix copy = this.copy();
			copy.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
				@Override
				public double visit(int row, int column, double value) {
					return FastMath.abs(value) < 1.0e-10 ? 0 : value;
				}
			});
			EigenDecompositionSymmetric eig = new EigenDecompositionSymmetric(copy, 1e-10, false);

			double[] eigenValues = eig.getEigenvalues();
			RealMatrix eigenVectors = eig.getVT();

			LOG.trace("eigenVectors " + eigenVectors);

			this.sHalf = MatrixUtils.createRealIdentityMatrix(this.getRowDimension());

			for (int i = 0; i < this.getRowDimension(); i++) {
				sHalf.setEntry(i, i, (sHalf.getEntry(i, i) / FastMath.sqrt(eigenValues[i])));
			}

			this.sHalf = eigenVectors.transpose().multiply(sHalf).multiply(eigenVectors);
		}

		return this.sHalf;
	}

	/** The SCF method used to compute overlap derivatives. */
	private SCFMethod scfMethod;

	/** The atom index with respect to which overlap derivatives are computed. */
	private int atomIndex;

	/**
	 * Compute Overlap partial derivative for an atom index.
	 * 
	 * @param atomIndex
	 *            the atom index with respect to which the derivative are to be
	 *            evaluated
	 * @param scfMethod
	 *            the reference to the SCFMethod
	 * @return three element array of Overlap elements representing partial
	 *         derivatives with respect to x, y and z of atom position
	 */
	public List<Overlap> computeDerivative(int atomIndex, SCFMethod scfMethod) {
		this.scfMethod = scfMethod;
		this.atomIndex = atomIndex;

		ArrayList<Overlap> dOverlap = new ArrayList<>(3);

		int noOfBasisFunctions = this.getRowDimension();
		Overlap dOverlapDx = new Overlap(noOfBasisFunctions);
		Overlap dOverlapDy = new Overlap(noOfBasisFunctions);
		Overlap dOverlapDz = new Overlap(noOfBasisFunctions);

		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < noOfBasisFunctions; j++) {
				Vector3D dOvrEle = computeOverlapDerElement(i, j);

				dOverlapDx.setEntry(i, j, dOvrEle.getX());
				dOverlapDy.setEntry(i, j, dOvrEle.getY());
				dOverlapDz.setEntry(i, j, dOvrEle.getZ());
			}
		}

		dOverlap.add(dOverlapDx);
		dOverlap.add(dOverlapDy);
		dOverlap.add(dOverlapDz);

		return dOverlap;
	}

	/**
	 * Compute one element of the overlap derivative matrix.
	 *
	 * @param i the row index (bra basis function)
	 * @param j the column index (ket basis function)
	 * @return the (i,j) element of the overlap derivative as a Vector3D (x, y, z components)
	 */
	private Vector3D computeOverlapDerElement(int i, int j) {
		BasisSetLibrary bsl = scfMethod.getOneEI().getBasisSetLibrary();
		ContractedGaussian cgi = bsl.getBasisFunctions().get(i);
		ContractedGaussian cgj = bsl.getBasisFunctions().get(j);

		return cgi.overlapDerivative(atomIndex, cgj);
	}

	/**
	 * Returns a string representation of this Overlap matrix, showing all elements.
	 *
	 * @return a formatted matrix string
	 */
	@Override
	public String toString() {
		return MathUtil.matrixToString(this);
	}
}
