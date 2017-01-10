package name.mjw.jquante.math.qm;

import java.util.ArrayList;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.Vector3D;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;

/**
 * Represents the HCore matrix.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class HCore extends Matrix {

	protected int atomIndex;
	private SCFMethod scfMethod;

	/**
	 * Creates a new instance of NxM Matrix
	 * 
	 * @param n
	 *            the first dimension
	 * @param m
	 *            the second dimension
	 */
	public HCore(int n, int m) {
		super(n, m);
	}

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public HCore(int n) {
		super(n, n);
	}

	/**
	 * Creates a new instance of Matrix, based on already allocated 2D array
	 * 
	 * @param a
	 *            the 2D array
	 */
	public HCore(double[][] a) {
		super(a);
	}

	/**
	 * Compute HCore partial derivative for an atom index.
	 * 
	 * @param atomIndex
	 *            the atom index with respect to which the derivative are to be
	 *            evaluated
	 * @param scfMethod
	 *            the reference to the SCFMethod
	 * @return three element array of HCore elements representing partial
	 *         derivatives with respect to x, y and z of atom position
	 */
	public ArrayList<HCore> computeDerivative(int atomIndex, SCFMethod scfMethod) {
		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;

		ArrayList<HCore> dHCore = new ArrayList<>(3);

		int noOfBasisFunctions = this.getRowCount();
		HCore dHCoreDx = new HCore(noOfBasisFunctions);
		HCore dHCoreDy = new HCore(noOfBasisFunctions);
		HCore dHCoreDz = new HCore(noOfBasisFunctions);

		double[][] hdx = dHCoreDx.getMatrix();
		double[][] hdy = dHCoreDy.getMatrix();
		double[][] hdz = dHCoreDz.getMatrix();

		int i;
		int j;
		for (i = 0; i < noOfBasisFunctions; i++) {
			for (j = 0; j < noOfBasisFunctions; j++) {
				Vector3D dHCoreEle = computeHCoreDerElement(atomIndex, i, j);

				hdx[i][j] = dHCoreEle.getI();
				hdy[i][j] = dHCoreEle.getJ();
				hdz[i][j] = dHCoreEle.getK();
			}
		}

		dHCore.add(dHCoreDx);
		dHCore.add(dHCoreDy);
		dHCore.add(dHCoreDz);

		return dHCore;
	}

	/**
	 * Compute one of the HCore derivative elements, with respect to an
	 * atomIndex
	 */
	private Vector3D computeHCoreDerElement(int atomIndex, int i, int j) {
		BasisFunctions bfs = scfMethod.getOneEI().getBasisFunctions();
		ContractedGaussian cgi = bfs.getBasisFunctions().get(i);
		ContractedGaussian cgj = bfs.getBasisFunctions().get(j);

		return cgi.kineticDerivative(atomIndex, cgj).add(
				cgi.nuclearAttractionDerivative(scfMethod.getMolecule(),
						atomIndex, cgj));

	}
}
