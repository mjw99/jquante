package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.test.Fixtures;

/**
 * Coverage for two-electron integral evaluation and Cauchy-Schwarz screening
 * with higher angular momentum basis functions. The rest of the suite uses
 * STO-3G only (s and p), so d and f shells - and therefore the screening and
 * packed-index code paths at high angular momentum - were previously untested.
 * d functions are taken from 6-31G** and f functions from cc-pVTZ.
 */
class HigherAngularMomentumIntegralsTest {

	private static int firstIndexWithTotalL(List<ContractedGaussian> bfs, int l) {
		for (int i = 0; i < bfs.size(); i++) {
			if (bfs.get(i).getTotalAngularMomentum() == l) {
				return i;
			}
		}
		throw new IllegalStateException("no basis function with total angular momentum " + l);
	}

	/**
	 * Asserts the fundamental properties every two-electron integral must satisfy:
	 * it is finite, obeys the eight-fold permutational symmetry, and respects the
	 * Cauchy-Schwarz bound used by the screening.
	 */
	private static void assertQuartetWellBehaved(TwoElectronIntegrals tei, int i, int j, int k, int l) {
		double v = tei.compute2E(i, j, k, l);
		assertTrue(Double.isFinite(v), "integral (" + i + j + "|" + k + l + ") must be finite");

		// The eight equivalent permutations of (ij|kl) must all agree.
		assertEquals(v, tei.compute2E(j, i, k, l), 1e-10);
		assertEquals(v, tei.compute2E(i, j, l, k), 1e-10);
		assertEquals(v, tei.compute2E(j, i, l, k), 1e-10);
		assertEquals(v, tei.compute2E(k, l, i, j), 1e-10);
		assertEquals(v, tei.compute2E(l, k, i, j), 1e-10);
		assertEquals(v, tei.compute2E(k, l, j, i), 1e-10);
		assertEquals(v, tei.compute2E(l, k, j, i), 1e-10);

		// Cauchy-Schwarz: |(ij|kl)| <= sqrt((ij|ij)) * sqrt((kl|kl)).
		double qij = Math.sqrt(Math.abs(tei.compute2E(i, j, i, j)));
		double qkl = Math.sqrt(Math.abs(tei.compute2E(k, l, k, l)));
		assertTrue(Math.abs(v) <= qij * qkl + 1e-9,
				"Schwarz bound violated for (" + i + j + "|" + k + l + "): |v|=" + Math.abs(v) + " bound=" + qij * qkl);
	}

	@Test
	void dFunctionIntegralsAreWellBehaved() throws Exception {
		BasisSetLibrary bsl = new BasisSetLibrary(Fixtures.getWater(), "6-31gss");
		List<ContractedGaussian> bfs = bsl.getBasisFunctions();
		int s = firstIndexWithTotalL(bfs, 0);
		int p = firstIndexWithTotalL(bfs, 1);
		int d = firstIndexWithTotalL(bfs, 2);
		TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true); // on-the-fly

		assertQuartetWellBehaved(tei, d, d, d, d);
		assertQuartetWellBehaved(tei, d, s, d, s);
		assertQuartetWellBehaved(tei, d, p, d, p);
		assertQuartetWellBehaved(tei, d, d, p, s);
	}

	@Test
	void fFunctionIntegralsAreWellBehaved() throws Exception {
		BasisSetLibrary bsl = new BasisSetLibrary(Fixtures.getHydrogenFluoride(), "cc-pvtz");
		List<ContractedGaussian> bfs = bsl.getBasisFunctions();
		int s = firstIndexWithTotalL(bfs, 0);
		int p = firstIndexWithTotalL(bfs, 1);
		int d = firstIndexWithTotalL(bfs, 2);
		int f = firstIndexWithTotalL(bfs, 3);
		TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true); // on-the-fly

		assertQuartetWellBehaved(tei, f, f, f, f);
		assertQuartetWellBehaved(tei, f, s, f, s);
		assertQuartetWellBehaved(tei, f, d, f, d);
		assertQuartetWellBehaved(tei, f, f, d, d);
		assertQuartetWellBehaved(tei, f, p, d, s);
	}

	@Test
	void screenedInCoreMatchesOnTheFlyWithDFunctions() throws Exception {
		// Drives the screened in-core compute2E path and packed indexing at d level
		// (water/6-31G** has 6 d functions); every stored integral must match the
		// unscreened on-the-fly value within the screening threshold.
		BasisSetLibrary bsl = new BasisSetLibrary(Fixtures.getWater(), "6-31gss");
		TwoElectronIntegrals inCore = new TwoElectronIntegrals(bsl); // screened, stored
		TwoElectronIntegrals onTheFly = new TwoElectronIntegrals(bsl, true); // unscreened singles
		double[] ints = inCore.getTwoEIntegrals();
		int n = bsl.getBasisFunctions().size();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++) {
				int ij = i * (i + 1) / 2 + j;
				for (int k = 0; k < n; k++) {
					for (int l = 0; l <= k; l++) {
						int kl = k * (k + 1) / 2 + l;
						if (ij >= kl) {
							double stored = ints[IntegralsUtil.ijkl2intindex(i, j, k, l)];
							assertEquals(onTheFly.compute2E(i, j, k, l), stored, 1e-9);
						}
					}
				}
			}
		}
	}
}
