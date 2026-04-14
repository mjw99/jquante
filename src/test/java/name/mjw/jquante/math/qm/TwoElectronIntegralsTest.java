package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

class TwoElectronIntegralsTest {

    private static final double DIFF = 1e-6;
    static BasisSetLibrary bsl;
    static Molecule hydrogen;

    @BeforeAll
    static void setUp() throws Exception {
        hydrogen = Fixtures.getHydrogen();
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
    void shellPairMethodMatchesConventionalForH2STO3G() {
        TwoElectronIntegrals conventional = new TwoElectronIntegrals(bsl);
        TwoElectronIntegrals shellPair = new TwoElectronIntegrals(bsl, hydrogen, false);

        assertFalse(shellPair.isOnTheFly());
        assertNotNull(shellPair.getTwoEIntegrals());
        assertEquals(conventional.getTwoEIntegrals().length, shellPair.getTwoEIntegrals().length);
        assertArrayEquals(conventional.getTwoEIntegrals(), shellPair.getTwoEIntegrals(), DIFF,
                "Shell-pair integrals must match conventional integrals for H2/STO-3G");
    }

    @Test
    void shellPairMethodMatchesConventionalForWaterSTO3G() throws Exception {
        Molecule water = Fixtures.getWater();
        BasisSetLibrary waterBsl = new BasisSetLibrary(water, "sto-3g");

        TwoElectronIntegrals conventional = new TwoElectronIntegrals(waterBsl);
        TwoElectronIntegrals shellPair = new TwoElectronIntegrals(waterBsl, water, false);

        assertFalse(shellPair.isOnTheFly());
        assertEquals(conventional.getTwoEIntegrals().length, shellPair.getTwoEIntegrals().length);
        assertArrayEquals(conventional.getTwoEIntegrals(), shellPair.getTwoEIntegrals(), DIFF,
                "Shell-pair integrals must match conventional integrals for H2O/STO-3G");
    }

    @Test
    void shellPairOnTheFlyDoesNotStoreIntegrals() {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, hydrogen, true);
        assertTrue(tei.isOnTheFly());
        assertNull(tei.getTwoEIntegrals());
    }

}
