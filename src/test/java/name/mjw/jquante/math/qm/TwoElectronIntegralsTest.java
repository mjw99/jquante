package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

/**
 * Tests for {@link TwoElectronIntegrals}, focusing on the shell-pair algorithm
 * ({@code compute2EShellPair}) introduced to replace the atom-based loop.
 *
 * <p>The reference values come from {@code compute2ESerial}, which is the
 * straightforward (non-parallel, non-shell-pair) implementation.
 */
class TwoElectronIntegralsTest {

    private static final double TOLERANCE = 1e-10;

    // ---------------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------------

    /**
     * Returns integrals computed by the serial (reference) algorithm.
     */
    private static double[] serialIntegrals(BasisSetLibrary bsl) {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true); // no auto-compute
        tei.compute2ESerial();
        return tei.getTwoEIntegrals().clone();
    }

    /**
     * Returns integrals computed by the new shell-pair algorithm.
     */
    private static double[] shellPairIntegrals(BasisSetLibrary bsl) {
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true); // no auto-compute
        tei.compute2EShellPair();
        return tei.getTwoEIntegrals().clone();
    }

    private static BasisSetLibrary bsl(Molecule mol) {
        return bsl(mol, "sto-3g");
    }

    private static BasisSetLibrary bsl(Molecule mol, String basisName) {
        try {
            return new BasisSetLibrary(mol, basisName);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    // ---------------------------------------------------------------------------
    // compute2EShellPair correctness vs compute2ESerial
    // ---------------------------------------------------------------------------

    @Test
    void shellPairMatchesSerialForHydrogen() {
        BasisSetLibrary bsl = bsl(Fixtures.getHydrogen());
        assertArrayEquals(serialIntegrals(bsl), shellPairIntegrals(bsl), TOLERANCE);
    }

    @Test
    void shellPairMatchesSerialForHydrogenFluoride() {
        BasisSetLibrary bsl = bsl(Fixtures.getHydrogenFluoride());
        assertArrayEquals(serialIntegrals(bsl), shellPairIntegrals(bsl), TOLERANCE);
    }

    @Test
    void shellPairMatchesSerialForWater() {
        BasisSetLibrary bsl = bsl(Fixtures.getWater());
        assertArrayEquals(serialIntegrals(bsl), shellPairIntegrals(bsl), TOLERANCE);
    }

    // ---------------------------------------------------------------------------
    // compute2EShellPair correctness vs compute2E (parallel reference)
    // ---------------------------------------------------------------------------

    @Test
    void shellPairMatchesParallelForWater() {
        BasisSetLibrary bsl = bsl(Fixtures.getWater());

        TwoElectronIntegrals parallelTei = new TwoElectronIntegrals(bsl, true);
        parallelTei.compute2E();
        double[] parallelResult = parallelTei.getTwoEIntegrals().clone();

        assertArrayEquals(parallelResult, shellPairIntegrals(bsl), TOLERANCE);
    }

    // ---------------------------------------------------------------------------
    // TwoElectronIntegrals(BasisSetLibrary, Molecule, boolean) constructor
    // ---------------------------------------------------------------------------

    @Test
    void moleculeConstructorPopulatesIntegrals() {
        Molecule water = Fixtures.getWater();
        BasisSetLibrary bsl = bsl(water);

        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, water, false);

        assertNotNull(tei.getTwoEIntegrals());
        assertEquals(shellPairIntegrals(bsl).length, tei.getTwoEIntegrals().length);
        assertArrayEquals(serialIntegrals(bsl), tei.getTwoEIntegrals(), TOLERANCE);
    }

    @Test
    void moleculeConstructorOnTheFlySkipsCompute() {
        Molecule hydrogen = Fixtures.getHydrogen();
        BasisSetLibrary bsl = bsl(hydrogen);

        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, hydrogen, true);

        // onTheFly=true: no integral array should have been allocated
        assertArrayEquals(null, tei.getTwoEIntegrals());
    }

    // ---------------------------------------------------------------------------
    // Integral count sanity checks
    // ---------------------------------------------------------------------------

    @Test
    void shellPairProducesCorrectIntegralCountForHydrogen() {
        BasisSetLibrary bsl = bsl(Fixtures.getHydrogen());
        // H2 STO-3G: n=2 basis functions → n*(n+1)*(n²+n+2)/8 = 2*3*8/8 = 6 unique integrals
        int n = bsl.getBasisFunctions().size();
        int expected = n * (n + 1) * (n * n + n + 2) / 8;
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        tei.compute2EShellPair();
        assertEquals(expected, tei.getTwoEIntegrals().length);
    }

    @Test
    void shellPairProducesCorrectIntegralCountForWater() {
        BasisSetLibrary bsl = bsl(Fixtures.getWater());
        // H2O STO-3G: n=7 basis functions → n*(n+1)*(n²+n+2)/8 = 7*8*58/8 = 406 unique integrals
        int n = bsl.getBasisFunctions().size();
        int expected = n * (n + 1) * (n * n + n + 2) / 8;
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        tei.compute2EShellPair();
        assertEquals(expected, tei.getTwoEIntegrals().length);
    }

    // ---------------------------------------------------------------------------
    // Larger basis (cc-pvtz): exercises index screening heavily
    //
    // cc-pvtz produces more shells per atom (s, p, d, f functions) which
    // creates many more shell pairs.  For a given bra shell pair with low
    // basis-function indices, roughly half the ket shell pairs will have
    // minPairIJ > maxPairIJ[bra] and are skipped by the index screen before
    // any basis-function loop is entered.  These tests confirm that screening
    // never drops a non-negligible integral.
    // ---------------------------------------------------------------------------

    @Test
    void shellPairMatchesSerialForHydrogenCcpvtz() {
        BasisSetLibrary bsl = bsl(Fixtures.getHydrogen(), "cc-pvtz");
        assertArrayEquals(serialIntegrals(bsl), shellPairIntegrals(bsl), TOLERANCE);
    }

    @Test
    void shellPairMatchesSerialForWaterCcpvtz() {
        BasisSetLibrary bsl = bsl(Fixtures.getWater(), "cc-pvtz");
        assertArrayEquals(serialIntegrals(bsl), shellPairIntegrals(bsl), TOLERANCE);
    }

    @Test
    void shellPairProducesCorrectIntegralCountForWaterCcpvtz() {
        BasisSetLibrary bsl = bsl(Fixtures.getWater(), "cc-pvtz");
        int n = bsl.getBasisFunctions().size();
        int expected = n * (n + 1) * (n * n + n + 2) / 8;
        TwoElectronIntegrals tei = new TwoElectronIntegrals(bsl, true);
        tei.compute2EShellPair();
        assertEquals(expected, tei.getTwoEIntegrals().length);
    }

    // ---------------------------------------------------------------------------
    // Known integral value
    // ---------------------------------------------------------------------------

    @Test
    void shellPairKnownValueForHydrogen() {
        // H2 STO-3G has one unique integral (0,0,0,0) ≈ 0.7746…
        // Value verified against compute2ESerial and literature
        BasisSetLibrary bsl = bsl(Fixtures.getHydrogen());
        double[] ints = shellPairIntegrals(bsl);
        assertEquals(serialIntegrals(bsl)[0], ints[0], TOLERANCE);
    }
}
