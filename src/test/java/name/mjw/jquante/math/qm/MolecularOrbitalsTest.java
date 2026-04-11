package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

class MolecularOrbitalsTest {

    private static final double DIFF = 0.001;

    static MolecularOrbitals mosFromHCore;
    static OneElectronIntegrals e1;
    static BasisSetLibrary bsl;

    @BeforeAll
    static void setUp() throws Exception {
        Molecule hydrogen = Fixtures.getHydrogen();
        bsl = new BasisSetLibrary(hydrogen, "sto-3g");
        e1 = new OneElectronIntegrals(bsl, hydrogen);

        mosFromHCore = new MolecularOrbitals(2);
        mosFromHCore.compute(e1.getHCore(), e1.getOverlap());
    }

    @Test
    void dimensionMatchesBasisSetSize() {
        assertEquals(2, mosFromHCore.getRowDimension());
        assertEquals(2, mosFromHCore.getColumnDimension());
    }

    @Test
    void getCoefficientsReturnsSquareMatrix() {
        double[][] coeffs = mosFromHCore.getCoefficients();
        assertNotNull(coeffs);
        assertEquals(2, coeffs.length);
        assertEquals(2, coeffs[0].length);
    }

    @Test
    void getOrbitalEnergiesReturnsCorrectCount() {
        double[] ev = mosFromHCore.getOrbitalEnergies();
        assertNotNull(ev);
        assertEquals(2, ev.length);
    }

    @Test
    void orbitalEnergiesFromHCoreAreOrdered() {
        // Lowest eigenvalue (most stabilizing) should come first
        double[] ev = mosFromHCore.getOrbitalEnergies();
        assertTrue(ev[0] <= ev[1],
                "First orbital energy should be <= second: ev[0]=" + ev[0] + " ev[1]=" + ev[1]);
    }

    @Test
    void computeFromConvergedFockGivesKnownOrbitalEnergies() throws Exception {
        // Using converged H2/STO-3G density to build the Fock, then diagonalize.
        // Expected orbital energies match the SCF-converged values.
        Molecule hydrogen = Fixtures.getHydrogen();
        BasisSetLibrary bf = new BasisSetLibrary(hydrogen, "sto-3g");
        OneElectronIntegrals oneEI = new OneElectronIntegrals(bf, hydrogen);
        TwoElectronIntegrals twoEI = new TwoElectronIntegrals(bf);

        // Converged density for H2/STO-3G (from GMatrixTest reference values)
        Density density = new Density(new double[][]{
                {0.3012278366, 0.3012278366},
                {0.3012278366, 0.3012278366}});

        GMatrix gMatrix = new GMatrix(2);
        gMatrix.compute(SCFType.HARTREE_FOCK, twoEI, density);

        Fock fock = new Fock(2);
        fock.compute(oneEI.getHCore(), gMatrix);

        MolecularOrbitals mos = new MolecularOrbitals(2);
        mos.compute(fock, oneEI.getOverlap());

        double[] ev = mos.getOrbitalEnergies();
        assertEquals(-0.578554, ev[0], DIFF);
        assertEquals(0.671144, ev[1], DIFF);
    }
}
