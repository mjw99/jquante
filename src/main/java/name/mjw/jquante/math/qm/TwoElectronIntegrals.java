package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PrimitiveGaussian;
import name.mjw.jquante.math.qm.basis.Shell;
import name.mjw.jquante.math.qm.integral.Integrals;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.molecule.Molecule;

import net.jafama.FastMath;

/**
 * The 2E integral driver.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class TwoElectronIntegrals {
	private static final Logger LOG = LogManager.getLogger(TwoElectronIntegrals.class);

	private BasisSetLibrary basisSetLibrary;

	/**
	 * Holds value of property twoEIntegrals.
	 */
	private double[] twoEIntegrals;

	/** The index of the atom with respect to which 2E derivatives are computed. */
	private int atomIndex;

	/** The SCF method instance used during derivative evaluation. */
	protected SCFMethod scfMethod;

	/** List of partial derivative arrays [dx, dy, dz] of the two-electron integrals. */
	protected ArrayList<double[]> twoEDer;

	/** If true, integrals are computed on-the-fly rather than stored in memory. */
	protected boolean onTheFly;

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisSetLibrary the basis functions to be used
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary) {
		this(basisSetLibrary, false);
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 *
	 * @param basisSetLibrary the basis functions to be used
	 * @param onTheFly        if true, the 2E integrals are not calculated and
	 *                        stored, they must be calculated individually by
	 *                        calling compute2E(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary, boolean onTheFly) {
		this.basisSetLibrary = basisSetLibrary;

		this.onTheFly = onTheFly;

		if (!onTheFly) {
			try {
				compute2EShellPair();
			} catch (OutOfMemoryError e) {
				LOG.error("No memory for in-core integral evaluation. Switching to direct integral evaluation.");
				this.onTheFly = true;
			}
		}
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals. The {@code molecule}
	 * parameter is accepted for API compatibility but is no longer required by
	 * the shell-pair algorithm; this constructor delegates to
	 * {@link #TwoElectronIntegrals(BasisSetLibrary, boolean)}.
	 *
	 * @param basisSetLibrary Basis functions for the molecule.
	 * @param molecule        Unused; kept for backwards compatibility.
	 * @param onTheFly        if true, the 2E integrals are not calculated and
	 *                        stored, they must be calculated individually by
	 *                        calling compute2E(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary, Molecule molecule, boolean onTheFly) {
		this(basisSetLibrary, onTheFly);
	}

	/**
	 * Compute the 2E integrals, and store it in a single 1D array, in the form
	 * [ijkl].
	 *
	 */
	protected void compute2ESerial() {
		LOG.debug("compute2ESerial() called");
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		// allocate required memory
		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		LOG.debug("noOfIntegrals is {}", noOfIntegrals);

		twoEIntegrals = new double[noOfIntegrals];

		// we only need i <= j, k <= l, and ij >= kl
		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < i + 1; j++) {
				int ij = i * (i + 1) / 2 + j;

				for (int k = 0; k < noOfBasisFunctions; k++) {
					for (int l = 0; l < (k + 1); l++) {
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
							// record the 2E integrals
							twoEIntegrals[ijkl] = Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
						}
					}
				}
			}
		}

	}

	/**
	 * Compute the 2E integrals, and store it in a single 1D array, in the form
	 * [ijkl].
	 * 
	 * This method has been modified to take advantage of multi core systems where
	 * available.
	 */
	protected void compute2E() {
		LOG.debug("compute2E() called");
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		// allocate required memory
		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		LOG.debug("noOfIntegrals is {}", noOfIntegrals);

		twoEIntegrals = new double[noOfIntegrals];

		// we only need i <= j, k <= l, and ij >= kl
		IntStream.range(0, noOfBasisFunctions).parallel().forEach(i -> {
			IntStream.range(0, i + 1).parallel().forEach(j -> {
				int ij = i * (i + 1) / 2 + j;
				IntStream.range(0, noOfBasisFunctions).parallel().forEach(k -> {
					IntStream.range(0, k + 1).parallel().forEach(l -> {
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
							// record the 2E integrals
							twoEIntegrals[ijkl] = Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
						}
					});
				});
			});
		});

	}

	/** Compute a single 2E derivative element */
	private Vector3D compute2EDerivativeElement(ContractedGaussian bfi, ContractedGaussian bfj, ContractedGaussian bfk,
			ContractedGaussian bfl) {

		Vector3D twoEDerEle = new Vector3D(0, 0, 0);
		int[] paramIdx;
		Power currentPower;
		Vector3D currentOrigin;
		double currentAlpha;

		if (bfi.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 4, 1, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				currentOrigin = iPG.origin();
				currentPower = iPG.powers();
				currentAlpha = iPG.exponent();
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						}
					}
				}
			}
		}

		if (bfj.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 4, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					currentOrigin = jPG.origin();
					currentPower = jPG.powers();
					currentAlpha = jPG.exponent();
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		if (bfk.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 1, 4, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						currentOrigin = kPG.origin();
						currentPower = kPG.powers();
						currentAlpha = kPG.exponent();
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		if (bfl.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 1, 2, 4 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							currentOrigin = lPG.origin();
							currentPower = lPG.powers();
							currentAlpha = lPG.exponent();

							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		return twoEDerEle;
	}

	/** helper for computing 2E derivative */
	private void compute2EDerivativeElementHelper(PrimitiveGaussian iPG, PrimitiveGaussian jPG, PrimitiveGaussian kPG,
			PrimitiveGaussian lPG, Vector3D currentOrigin, Power currentPower, double currentAlpha, int[] paramIdx,
			Vector3D derEle) {

		int l = currentPower.l();
		int m = currentPower.m();
		int n = currentPower.n();

		double coeff = iPG.coefficient() * jPG.coefficient() * kPG.coefficient() * lPG.coefficient();

		PrimitiveGaussian xPG = new PrimitiveGaussian(currentOrigin, new Power(l + 1, m, n), currentAlpha, coeff);

		PrimitiveGaussian[] pgs = new PrimitiveGaussian[] { iPG, jPG, kPG, lPG, xPG };

		double terma = FastMath.sqrt(currentAlpha * (2.0 * l + 1.0)) * coeff
				* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		double termbx = 0.0;
		double termby = 0.0;
		double termbz = 0.0;

		if (l > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(xPG.powers().l() - 1, xPG.powers().m(), xPG.powers().n()), coeff, currentAlpha);
			termbx = -2.0 * l * FastMath.sqrt(currentAlpha / (2. * l - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		if (m > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(xPG.powers().l(), xPG.powers().m() - 1, xPG.powers().n()), coeff, currentAlpha);
			termby = -2.0 * m * FastMath.sqrt(currentAlpha / (2. * m - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		if (n > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(xPG.powers().l(), xPG.powers().m(), xPG.powers().n() - 1), coeff, currentAlpha);
			termbz = -2.0 * n * FastMath.sqrt(currentAlpha / (2. * n - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		derEle = new Vector3D(derEle.getX() + terma + termbx, derEle.getY() + terma + termby,
				derEle.getZ() + terma + termbz);
	}

	/**
	 * Compute the 2E integrals using a true shell-pair based method, and store
	 * them in a single 1D array in the form [ijkl].
	 *
	 * <p>Shell pairs (IJ) and (KL) are obtained from
	 * {@link BasisSetLibrary#getUniqueShellPairs()}. Two screening layers are
	 * applied before entering the basis-function loops:
	 * <ol>
	 *   <li><b>Index screening</b>: skip a quartet when
	 *       {@code maxPairIJ[bra] < minPairIJ[ket]}, i.e. every possible ij
	 *       in the bra is less than every possible kl in the ket, so
	 *       {@code ij >= kl} can never hold.</li>
	 *   <li><b>Schwarz screening</b>: skip a quartet when
	 *       {@code Q_IJ * Q_KL < threshold}, where
	 *       {@code Q_IJ = max_{i∈I,j∈J} sqrt(|(ij|ij)|)}.</li>
	 * </ol>
	 * The outer loop over bra shell pairs is parallelised.
	 */
	protected void compute2EShellPair() {
		LOG.debug("compute2EShellPair() called");
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		twoEIntegrals = new double[noOfIntegrals];

		final List<List<Shell>> shellPairs = basisSetLibrary.getUniqueShellPairs();
		final int nShellPairs = shellPairs.size();

		// Pre-compute per-shell-pair data into flat arrays to avoid repeated
		// list lookups and helper-method calls inside the hot inner loop.
		final int[] pHighFirst = new int[nShellPairs];
		final int[] pHighLast  = new int[nShellPairs];
		final int[] pLowFirst  = new int[nShellPairs];
		final int[] pLowLast   = new int[nShellPairs];
		// minPairIJ / maxPairIJ: the minimum and maximum ij = i*(i+1)/2+j
		// achievable for valid (i,j) pairs within each shell pair.
		// For a cross-pair, j < i always, so min = iFirst*(iFirst+1)/2 + jFirst
		// and max = iLast*(iLast+1)/2 + jLast.
		// The same formula holds for self-pairs (iFirst==jFirst, iLast==jLast).
		final int[] minPairIJ   = new int[nShellPairs];
		final int[] maxPairIJ   = new int[nShellPairs];
		final double[] schwarzValues = new double[nShellPairs];

		for (int p = 0; p < nShellPairs; p++) {
			final Shell high = higherIndexedShell(shellPairs.get(p));
			final Shell low  = lowerIndexedShell(shellPairs.get(p));
			final int iFirst = high.getFirstBasisFunctionIndex();
			final int iLast  = high.getLastBasisFunctionIndex();
			final int jFirst = low.getFirstBasisFunctionIndex();
			final int jLast  = low.getLastBasisFunctionIndex();
			pHighFirst[p] = iFirst;
			pHighLast[p]  = iLast;
			pLowFirst[p]  = jFirst;
			pLowLast[p]   = jLast;
			minPairIJ[p] = iFirst * (iFirst + 1) / 2 + jFirst;
			maxPairIJ[p] = iLast  * (iLast  + 1) / 2 + jLast;

			// Schwarz bound: Q_p = max_{i,j in pair, j<=i} sqrt(|(ij|ij)|)
			double maxVal = 0.0;
			for (int i = iFirst; i <= iLast; i++) {
				for (int j = jFirst; j <= jLast; j++) {
					if (j > i) continue;
					double val = Math.abs(Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(i), bfs.get(j)));
					maxVal = Math.max(maxVal, Math.sqrt(val));
				}
			}
			schwarzValues[p] = maxVal;
		}

		final double SCHWARZ_THRESHOLD = 1.0e-10;

		// Parallelize over ALL (ijShell, klShell) pairs as a flat index.
		//
		// This gives P² parallel tasks (vs P tasks with a sequential inner loop),
		// matching the parallelism granularity of compute2E() while retaining the
		// Schwarz and index-screening benefits of the shell-pair approach.
		//
		// Correctness: each unique integral (i,j,k,l) is written by exactly one task.
		// Proof: BF indices uniquely determine their shell pair; if ijShell ≠ klShell
		// the BF sets are disjoint, so ij ≠ kl, and only the task where ij ≥ kl
		// passes the BF-level guard. The other direction fails the guard and writes
		// nothing. For ijShell == klShell the diagonal task handles both roles.
		IntStream.range(0, nShellPairs * nShellPairs).parallel().forEach(idx -> {
			final int ijShell = idx / nShellPairs;
			final int klShell = idx % nShellPairs;

			// Index screening: skip if every ij in bra < every kl in ket
			if (maxPairIJ[ijShell] < minPairIJ[klShell]) return;

			// Schwarz screening: skip if the integral upper bound is negligible
			if (schwarzValues[ijShell] * schwarzValues[klShell] < SCHWARZ_THRESHOLD) return;

			final int iFirst = pHighFirst[ijShell];
			final int iLast  = pHighLast[ijShell];
			final int jFirst = pLowFirst[ijShell];
			final int jLast  = pLowLast[ijShell];
			final int kFirst = pHighFirst[klShell];
			final int kLast  = pHighLast[klShell];
			final int lFirst = pLowFirst[klShell];
			final int lLast  = pLowLast[klShell];

			for (int i = iFirst; i <= iLast; i++) {
				for (int j = jFirst; j <= jLast; j++) {
					if (j > i) continue;
					final int ij = i * (i + 1) / 2 + j;

					for (int k = kFirst; k <= kLast; k++) {
						for (int l = lFirst; l <= lLast; l++) {
							if (l > k) continue;
							final int kl = k * (k + 1) / 2 + l;

							if (ij >= kl) {
								twoEIntegrals[IntegralsUtil.ijkl2intindex(i, j, k, l)] =
										Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
							}
						}
					}
				}
			}
		});
	}

	/**
	 * Returns the shell in the pair with the higher (or equal) first basis function
	 * index, which plays the role of the "i" index in the upper-triangular loop.
	 */
	private static Shell higherIndexedShell(List<Shell> pair) {
		return pair.get(0).getFirstBasisFunctionIndex() >= pair.get(1).getFirstBasisFunctionIndex()
				? pair.get(0) : pair.get(1);
	}

	/**
	 * Returns the shell in the pair with the lower (or equal) first basis function
	 * index, which plays the role of the "j" index in the upper-triangular loop.
	 */
	private static Shell lowerIndexedShell(List<Shell> pair) {
		return pair.get(0).getFirstBasisFunctionIndex() <= pair.get(1).getFirstBasisFunctionIndex()
				? pair.get(0) : pair.get(1);
	}

	/**
	 * 
	 * Compute an integral centered at &lt;ij&#x7C;kl&gt;
	 * 
	 * @param i Index of contracted Gaussian function i.
	 * @param j Index of contracted Gaussian function j.
	 * @param k Index of contracted Gaussian function k.
	 * @param l Index of contracted Gaussian function l.
	 * @return the value of two integral electron
	 */
	public double compute2E(int i, int j, int k, int l) {
		// if we really need to compute, go ahead!
		List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		return Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
	}

	/**
	 * Compute an integral centered at &lt;ij&#x7C;kl&gt;
	 * 
	 * @param cgi Contracted Gaussian function i.
	 * @param cgj Contracted Gaussian function j.
	 * @param cgk Contracted Gaussian function k.
	 * @param cgl Contracted Gaussian function l.
	 * @return the value of two integral electron
	 */
	public double compute2E(ContractedGaussian cgi, ContractedGaussian cgj, ContractedGaussian cgk,
			ContractedGaussian cgl) {
		return Integrals.coulomb(cgi, cgj, cgk, cgl);
	}

	/**
	 * Getter for property twoEIntegrals.
	 * 
	 * @return Value of property twoEIntegrals.
	 */
	public double[] getTwoEIntegrals() {
		return this.twoEIntegrals;
	}

	/**
	 * Setter for property twoEIntegrals.
	 * 
	 * @param twoEIntegrals New value of property twoEIntegrals.
	 */
	public void setTwoEIntegrals(double[] twoEIntegrals) {
		this.twoEIntegrals = twoEIntegrals;
	}

	/**
	 * Return the derivatives of the two electron integrals.
	 * 
	 * @param atomIndex the reference atom index
	 * @param scfMethod the reference SCF method
	 */
	public void compute2EDerivatives(int atomIndex, SCFMethod scfMethod) {
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();
		final int noOfBasisFunctions = bfs.size();

		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;

		twoEDer = new ArrayList<>();

		double[] dxTwoE = new double[twoEIntegrals.length];
		double[] dyTwoE = new double[twoEIntegrals.length];
		double[] dzTwoE = new double[twoEIntegrals.length];

		twoEDer.add(dxTwoE);
		twoEDer.add(dyTwoE);
		twoEDer.add(dzTwoE);

		// we only need i <= j, k <= l, and ij >= kl
		IntStream.range(0, noOfBasisFunctions).parallel().forEach(i -> {
			ContractedGaussian bfi = bfs.get(i);

			IntStream.range(0, i + 1).parallel().forEach(j -> {
				ContractedGaussian bfj = bfs.get(j);
				int ij = i * (i + 1) / 2 + j;
				IntStream.range(0, noOfBasisFunctions).parallel().forEach(k -> {
					ContractedGaussian bfk = bfs.get(k);
					IntStream.range(0, k + 1).parallel().forEach(l -> {
						ContractedGaussian bfl = bfs.get(l);
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);

							// record derivative of the 2E integrals
							Vector3D twoEDerEle = compute2EDerivativeElement(bfi, bfj, bfk, bfl);

							dxTwoE[ijkl] = twoEDerEle.getX();
							dyTwoE[ijkl] = twoEDerEle.getY();
							dzTwoE[ijkl] = twoEDerEle.getZ();

						}
					});
				});
			});
		});
	}

	/**
	 * The currently computed list of two electron derivatives
	 * 
	 * @return partial derivatives of 2E integrals, computed in previous call to
	 *         compute2EDerivatives()
	 */
	public ArrayList<double[]> getTwoEDer() {
		return twoEDer;
	}

	/**
	 * Get the value of onTheFly
	 * 
	 * @return the value of onTheFly
	 */
	public boolean isOnTheFly() {
		return onTheFly;
	}

	/**
	 * Set the value of onTheFly
	 * 
	 * @param onTheFly new value of onTheFly
	 */
	public void setOnTheFly(boolean onTheFly) {
		this.onTheFly = onTheFly;
	}

	/**
	 * Prints all two-electron integrals (i,j,k,l) and their values to standard output.
	 * Intended for debugging purposes.
	 */
	public void printCompleteIntegralList() {

		int noOfBasisFunctions = basisSetLibrary.getBasisFunctions().size();

		System.out.println("");
		System.out.println("Complete integral list");
		System.out.println("======================");
		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < noOfBasisFunctions; j++) {
				for (int k = 0; k < noOfBasisFunctions; k++) {
					for (int l = 0; l < noOfBasisFunctions; l++) {
						int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
						System.out.println(i + " " + j + " " + k + " " + l + "\t\t" + twoEIntegrals[ijkl]);

					}

				}

			}

		}
	}
}
