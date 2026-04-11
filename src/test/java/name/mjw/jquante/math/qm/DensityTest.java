package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class DensityTest {

    @Test
    void constructorNCreatesDimension() {
        Density d = new Density(3);
        assertEquals(3, d.getRowDimension());
        assertEquals(3, d.getColumnDimension());
    }

    @Test
    void constructorFromDataSetsEntries() {
        double[][] data = {{1.0, 0.5}, {0.5, 2.0}};
        Density d = new Density(data);
        assertEquals(1.0, d.getEntry(0, 0), 1e-10);
        assertEquals(0.5, d.getEntry(0, 1), 1e-10);
        assertEquals(0.5, d.getEntry(1, 0), 1e-10);
        assertEquals(2.0, d.getEntry(1, 1), 1e-10);
    }

    @Test
    void toStringReturnsNonEmpty() {
        Density d = new Density(new double[][]{{1.0, 0.0}, {0.0, 1.0}});
        String s = d.toString();
        assertNotNull(s);
        assertFalse(s.isEmpty());
    }

    @Test
    void computeFromMolecularOrbitalsWithOneOccupied() {
        // P = C_occ^T * C_occ where C_occ = first row of MO matrix
        // mos row 0 = [1.0, 0.0], noOfOccupiedMOs = 1
        // => P = [[1,0],[0,0]]
        MolecularOrbitals mos = new MolecularOrbitals(2);
        mos.setSubMatrix(new double[][]{{1.0, 0.0}, {0.0, 1.0}}, 0, 0);

        Density density = new Density(2);
        density.compute(null, false, null, 1, mos);

        assertEquals(1.0, density.getEntry(0, 0), 1e-10);
        assertEquals(0.0, density.getEntry(0, 1), 1e-10);
        assertEquals(0.0, density.getEntry(1, 0), 1e-10);
        assertEquals(0.0, density.getEntry(1, 1), 1e-10);
    }

    @Test
    void computeFromMolecularOrbitalsWithTwoOccupied() {
        // With noOfOccupiedMOs=2, both rows contribute
        // mos = identity, C_occ = identity, P = I^T * I = I
        MolecularOrbitals mos = new MolecularOrbitals(2);
        mos.setSubMatrix(new double[][]{{1.0, 0.0}, {0.0, 1.0}}, 0, 0);

        Density density = new Density(2);
        density.compute(null, false, null, 2, mos);

        assertEquals(1.0, density.getEntry(0, 0), 1e-10);
        assertEquals(1.0, density.getEntry(1, 1), 1e-10);
        assertEquals(0.0, density.getEntry(0, 1), 1e-10);
    }

    @Test
    void computeWithDensityGuesser() {
        DensityGuesser guesser = scf -> new Density(new double[][]{{3.0, 0.0}, {0.0, 3.0}});

        Density density = new Density(2);
        density.compute(null, true, guesser, 1, null);

        assertEquals(3.0, density.getEntry(0, 0), 1e-10);
        assertEquals(3.0, density.getEntry(1, 1), 1e-10);
        assertEquals(0.0, density.getEntry(0, 1), 1e-10);
    }

    @Test
    void guessInitialDmFalseWithNullGuesserUsesOrbitals() {
        MolecularOrbitals mos = new MolecularOrbitals(2);
        mos.setSubMatrix(new double[][]{{0.5, 0.5}, {0.5, 0.5}}, 0, 0);

        Density density = new Density(2);
        // guessInitialDM=true but guesser=null: condition is false, falls through to MO path
        density.compute(null, true, null, 1, mos);

        // C_occ = [0.5, 0.5], P = [[0.25, 0.25], [0.25, 0.25]]
        assertEquals(0.25, density.getEntry(0, 0), 1e-10);
        assertEquals(0.25, density.getEntry(0, 1), 1e-10);
    }
}
