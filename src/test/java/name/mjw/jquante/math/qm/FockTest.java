package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class FockTest {

    @Test
    void constructorNCreatesDimension() {
        Fock f = new Fock(4);
        assertEquals(4, f.getRowDimension());
        assertEquals(4, f.getColumnDimension());
    }

    @Test
    void constructorFromDataSetsEntries() {
        double[][] data = {{1.0, 2.0}, {3.0, 4.0}};
        Fock f = new Fock(data);
        assertEquals(1.0, f.getEntry(0, 0), 1e-10);
        assertEquals(2.0, f.getEntry(0, 1), 1e-10);
        assertEquals(3.0, f.getEntry(1, 0), 1e-10);
        assertEquals(4.0, f.getEntry(1, 1), 1e-10);
    }

    @Test
    void computeIsSumOfHCoreAndGMatrix() {
        HCore hCore = new HCore(2);
        hCore.setEntry(0, 0, 1.0);
        hCore.setEntry(0, 1, 0.5);
        hCore.setEntry(1, 0, 0.5);
        hCore.setEntry(1, 1, 2.0);

        GMatrix gMatrix = new GMatrix(new double[][]{{0.1, 0.2}, {0.2, 0.3}});

        Fock fock = new Fock(2);
        fock.compute(hCore, gMatrix);

        assertEquals(1.1, fock.getEntry(0, 0), 1e-10);
        assertEquals(0.7, fock.getEntry(0, 1), 1e-10);
        assertEquals(0.7, fock.getEntry(1, 0), 1e-10);
        assertEquals(2.3, fock.getEntry(1, 1), 1e-10);
    }

    @Test
    void computeWithZeroGMatrixEqualsHCore() {
        HCore hCore = new HCore(2);
        hCore.setEntry(0, 0, -1.5);
        hCore.setEntry(0, 1, -0.3);
        hCore.setEntry(1, 0, -0.3);
        hCore.setEntry(1, 1, -0.8);

        GMatrix gMatrix = new GMatrix(2);  // all zeros

        Fock fock = new Fock(2);
        fock.compute(hCore, gMatrix);

        assertEquals(-1.5, fock.getEntry(0, 0), 1e-10);
        assertEquals(-0.3, fock.getEntry(0, 1), 1e-10);
        assertEquals(-0.8, fock.getEntry(1, 1), 1e-10);
    }

    @Test
    void toStringReturnsNonEmpty() {
        Fock f = new Fock(new double[][]{{1.0, 0.0}, {0.0, 1.0}});
        String s = f.toString();
        assertNotNull(s);
        assertFalse(s.isEmpty());
    }
}
