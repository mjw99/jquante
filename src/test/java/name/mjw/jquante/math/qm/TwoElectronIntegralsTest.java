package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

class TwoElectronIntegralsTest {

    private static final double DIFF = 1e-6;
    static BasisSetLibrary bsl;

    @BeforeAll
    static void setUp() throws Exception {
        Molecule hydrogen = Fixtures.getHydrogen();
        bsl = new BasisSetLibrary(hydrogen, "sto-3g");
    }

    @Test
    void defaultConstructorComputesIntegralsInCore() {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl);
        assertFalse(tei.isOnTheFly());
        assertNotNull(tei.getTwoEIntegrals());
    }

    @Test
    void onTheFlyConstructorDoesNotStoreIntegrals() {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        assertTrue(tei.isOnTheFly());
        assertNull(tei.getTwoEIntegrals());
    }

    @Test
    void inCoreIntegralsCountForH2STO3G() {
        // For n=2 basis functions: n*(n+1)*(n^2+n+2)/8 = 2*3*8/8 = 6
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl);
        assertEquals(6, tei.getTwoEIntegrals().length);
    }

    @Test
    void coulombIntegralOfSameBasisFunctionIsPositive() {
        // (phi_0 phi_0 | phi_0 phi_0) is a Coulomb self-repulsion integral: must be > 0
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        double val = tei.compute2E(0, 0, 0, 0);
        assertTrue(val > 0, "Coulomb integral (0,0,0,0) should be positive, got: " + val);
    }

    @Test
    void compute2EByIndexMatchesStoredInCore() {
        TwoElectronIntegrals teiInCore = new TwoElectronIntegrals(bsl);
        TwoElectronIntegrals teiOnTheFly = new TwoElectronIntegrals(bsl, true);

        // (0,0,0,0) integral: stored == computed on-the-fly
        double stored = teiInCore.getTwoEIntegrals()[0];  // index 0 = ijkl(0,0,0,0)
        double onTheFly = teiOnTheFly.compute2E(0, 0, 0, 0);
        assertEquals(stored, onTheFly, DIFF);
    }

    @Test
    void setGetTwoEIntegrals() {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        double[] integrals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        tei.setTwoEIntegrals(integrals);
        assertArrayEquals(integrals, tei.getTwoEIntegrals(), DIFF);
    }

    @Test
    void setGetOnTheFly() {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        assertTrue(tei.isOnTheFly());
        tei.setOnTheFly(false);
        assertFalse(tei.isOnTheFly());
        tei.setOnTheFly(true);
        assertTrue(tei.isOnTheFly());
    }

    @Test
    void allInCoreIntegralsAreFinite() {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl);
        for (double v : tei.getTwoEIntegrals()) {
            assertTrue(Double.isFinite(v), "All stored integrals should be finite, found: " + v);
        }
    }

    @Test
    void integralStorageSizeMatchesUniqueCountForInCoreBasis() {
        assertEquals(6, TwoElectronIntegrals.integralStorageSize(2));
        // 350 functions is near the in-core limit but still addressable by an int
        assertEquals((int) IntegralsUtil.numberOfUniqueIntegrals(350),
                TwoElectronIntegrals.integralStorageSize(350));
    }

    @Test
    void integralStorageSizeThrowsWhenCountExceedsIntRange() {
        // 362 functions is the first size whose integral count exceeds
        // Integer.MAX_VALUE, so it cannot be stored in core.
        assertThrows(ArithmeticException.class, () -> TwoElectronIntegrals.integralStorageSize(362));
        assertThrows(ArithmeticException.class, () -> TwoElectronIntegrals.integralStorageSize(1000));
    }

    @Test
    void schwarzBoundIsAnUpperBoundForEveryQuartet() throws Exception {
        // Cauchy-Schwarz: |(ij|kl)| <= sqrt((ij|ij)) * sqrt((kl|kl)). The screening
        // is only correct if this holds for every quartet of the basis.
        BasisSetLibrary water = new BasisSetLibrary(Fixtures.getWater(), "sto-3g");
        TwoElectronIntegrals tei = new TwoElectronIntegrals(water, true); // on-the-fly, unscreened singles
        int n = water.getBasisFunctions().size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double qij = Math.sqrt(Math.abs(tei.compute2E(i, j, i, j)));
                for (int k = 0; k < n; k++) {
                    for (int l = 0; l <= k; l++) {
                        double qkl = Math.sqrt(Math.abs(tei.compute2E(k, l, k, l)));
                        double val = Math.abs(tei.compute2E(i, j, k, l));
                        double bound = qij * qkl;
                        assertTrue(val <= bound + 1e-9,
                                "Schwarz violated for (" + i + j + "|" + k + l + "): |val|=" + val + " bound=" + bound);
                    }
                }
            }
        }
    }

    @Test
    void defaultScreeningPreservesIntegralsForCompactSystem() throws Exception {
        // For a compact molecule nothing of significance is screened, so every
        // stored integral must still match the unscreened value (well within the
        // 1e-10 default threshold).
        BasisSetLibrary water = new BasisSetLibrary(Fixtures.getWater(), "sto-3g");
        TwoElectronIntegrals screened = new TwoElectronIntegrals(water);           // in-core, default screening
        TwoElectronIntegrals reference = new TwoElectronIntegrals(water, true);    // on-the-fly, unscreened
        double[] ints = screened.getTwoEIntegrals();
        int n = water.getBasisFunctions().size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                int ij = i * (i + 1) / 2 + j;
                for (int k = 0; k < n; k++) {
                    for (int l = 0; l <= k; l++) {
                        int kl = k * (k + 1) / 2 + l;
                        if (ij >= kl) {
                            double stored = ints[IntegralsUtil.ijkl2intindex(i, j, k, l)];
                            assertEquals(reference.compute2E(i, j, k, l), stored, 1e-9);
                        }
                    }
                }
            }
        }
    }

    @Test
    void aggressiveThresholdScreensEveryIntegral() {
        // With a huge threshold no quartet can clear the bound, so all are skipped
        // and left as zero.
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true); // H2, no auto-compute
        tei.setSchwarzThreshold(1.0e6);
        tei.compute2E();
        for (double v : tei.getTwoEIntegrals()) {
            assertEquals(0.0, v, 0.0);
        }
    }

    @Test
    void schwarzPredicateSkipsNegligibleQuartetsForBenzene() throws Exception {
        // Benzene/STO-3G is spread out enough that many distant quartets are
        // negligible. Apply the exact criterion the compute2E* methods use,
        // Q[ij]*Q[kl] < threshold, to the real integral magnitudes and confirm it
        // classifies a strict, non-empty subset of canonical quartets as skippable.
        BasisSetLibrary benzene = new BasisSetLibrary(Fixtures.getBenzene(), "sto-3g");
        TwoElectronIntegrals tei = new TwoElectronIntegrals(benzene, true);
        int n = benzene.getBasisFunctions().size();
        int nPairs = n * (n + 1) / 2;

        // Schwarz bounds Q[ij] = sqrt(|(ij|ij)|) -- only O(n^2) integral evaluations.
        double[] q = new double[nPairs];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                q[i * (i + 1) / 2 + j] = Math.sqrt(Math.abs(tei.compute2E(i, j, i, j)));
            }
        }

        double threshold = 1.0e-8;
        long total = 0;
        long skipped = 0;
        for (int ij = 0; ij < nPairs; ij++) {
            for (int kl = 0; kl <= ij; kl++) { // canonical quartets with ij >= kl
                total++;
                if (q[ij] * q[kl] < threshold) {
                    skipped++;
                }
            }
        }

        assertTrue(skipped > 0, "screening should skip negligible quartets, skipped=" + skipped);
        assertTrue(skipped < total, "screening must retain the significant quartets");
    }
}
