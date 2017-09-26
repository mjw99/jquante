package name.mjw.jquante.math.qm;

import java.util.ArrayList;

import org.apache.log4j.Logger;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.Vector;
import name.mjw.jquante.math.la.GaussianElimination;
import name.mjw.jquante.math.la.exception.SingularMatrixException;

/**
 * DIIS (Direct Inversion of iterative subspaces) proposed by Peter Pualay for
 * HF convergence acceleration.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class DIISFockExtrapolator implements FockExtrapolator {

	/**
	 * Logger object
	 */
	private static final Logger LOG = Logger
			.getLogger(DIISFockExtrapolator.class);

	private ArrayList<Fock> fockMatrixList = new ArrayList<>();
	private ArrayList<Vector> errorMatrixList = new ArrayList<>();

	protected int diisStep = 0;

	private boolean isDiisStarted = false;

	private Fock oldFock = null;

	/**
	 * Get the next extrapolated fock matrix
	 * 
	 * @param currentFock
	 *            the current fock
	 * @param overlap
	 *            the overlap matrix
	 * @param density
	 *            the current density matrix
	 * @return extrapolated fock
	 */
	@Override
	public Fock next(Fock currentFock, Overlap overlap, Density density) {
		Fock newFock = new Fock(((Matrix) currentFock.clone()).getMatrix());
		double[][] newFockMat = newFock.getMatrix();

		Matrix fPS = currentFock.mul(density).mul(overlap);
		Matrix sPF = overlap.mul(density).mul(currentFock);

		Vector errorMatrix = new Vector(fPS.sub(sPF));
		double mxerr = errorMatrix.maxNorm();

		if (mxerr < errorThreshold && !isDiisStarted) {
			isDiisStarted = true;
		}

		if (!isDiisStarted) {
			if (oldFock == null) {
				oldFock = currentFock;

				return currentFock;
			} else {
				newFockMat = oldFock.mul(0.5).add(currentFock.mul(0.5))
						.getMatrix();
				oldFock = currentFock;

				return newFock;
			}
		}

		errorMatrixList.add(errorMatrix);
		fockMatrixList.add(currentFock);

		newFock.makeZero();

		int noOfIterations = errorMatrixList.size();
		int n1 = noOfIterations + 1;

		Matrix a = new Matrix(n1);
		Vector b = new Vector(n1);

		double[][] aMatrix = a.getMatrix();
		double[] bVector = b.getVector();

		// set up A x = B to be solved
		for (int i = 0; i < noOfIterations; i++) {
			for (int j = 0; j < noOfIterations; j++) {
				aMatrix[i][j] = errorMatrixList.get(i).dot(
						errorMatrixList.get(j));
			}
		}

		for (int i = 0; i < noOfIterations; i++) {
			aMatrix[noOfIterations][i] = aMatrix[i][noOfIterations] = -1.0;
			bVector[i] = 0.0;
		}

		aMatrix[noOfIterations][noOfIterations] = 0.0;
		bVector[noOfIterations] = -1.0;

		GaussianElimination gele = new GaussianElimination();
		gele.setMatrixA(a);
		gele.setVectorX(b);
		try {
			gele.findSolution();
			double[] solVec = gele.getVectorX().getVector();

			LOG.debug("Solution found: " + gele.getVectorX());

			for (int i = 0; i < noOfIterations; i++) {
				double[][] prevFockMat = fockMatrixList.get(i).getMatrix();
				for (int j = 0; j < newFockMat.length; j++) {
					for (int k = 0; k < newFockMat.length; k++) {
						newFockMat[j][k] += solVec[i] * prevFockMat[j][k];
					}
				}
			}

			oldFock = currentFock;
		} catch (SingularMatrixException ignored) {
			// no solution could be found, so return the current Fock as is
			LOG.debug("No solution: " + diisStep);
			diisStep++;
			return currentFock;
		}

		diisStep++;

		return newFock;
	}

	protected double errorThreshold;

	/**
	 * Get the value of errorThreshold
	 * 
	 * @return the value of errorThreshold
	 */
	public double getErrorThreshold() {
		return errorThreshold;
	}

	/**
	 * Set the value of errorThreshold
	 * 
	 * @param errorThreshold
	 *            new value of errorThreshold
	 */
	public void setErrorThreshold(double errorThreshold) {
		this.errorThreshold = errorThreshold;
	}
}
