package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

class OverlapTest {

    private static final double DIFF = 1e-8;
    static Overlap overlap;

    @BeforeAll
    static void setUp() throws Exception {
        Molecule hydrogen = Fixtures.getHydrogen();
        BasisSetLibrary bsl = new BasisSetLibrary(hydrogen, "sto-3g");
        OneElectronIntegrals e1 = new OneElectronIntegrals(bsl, hydrogen);
        overlap = e1.getOverlap();
    }

    @Test
    void diagonalElementsAreOne() {
        assertEquals(1.0, overlap.getEntry(0, 0), DIFF);
        assertEquals(1.0, overlap.getEntry(1, 1), DIFF);
    }

    @Test
    void isSymmetric() {
        assertEquals(overlap.getEntry(0, 1), overlap.getEntry(1, 0), DIFF);
    }

    @Test
    void offDiagonalIsPositiveForHydrogenSTO3G() {
        // Two H 1s functions at bond distance should overlap positively
        assertTrue(overlap.getEntry(0, 1) > 0,
                "Off-diagonal overlap for H2/STO-3G should be positive");
    }

    @Test
    void sHalfIsInverseSquareRoot() {
        // X = S^{-1/2} satisfies X * S * X = I
        RealMatrix x = overlap.getSHalf();
        RealMatrix product = x.multiply(overlap).multiply(x);
        int n = overlap.getRowDimension();
        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                assertEquals(identity.getEntry(i, j), product.getEntry(i, j), 1e-6);
            }
        }
    }

    @Test
    void sHalfIsCached() {
        // Calling getSHalf() twice should return the same object instance
        RealMatrix x1 = overlap.getSHalf();
        RealMatrix x2 = overlap.getSHalf();
        assertSame(x1, x2);
    }

    @Test
    void sHalfIsSquareWithCorrectDimension() {
        RealMatrix x = overlap.getSHalf();
        assertEquals(overlap.getRowDimension(), x.getRowDimension());
        assertEquals(overlap.getColumnDimension(), x.getColumnDimension());
    }

    @Test
    void toStringReturnsNonEmpty() {
        String s = overlap.toString();
        assertNotNull(s);
        assertFalse(s.isEmpty());
    }

    @Test
    void constructorCreatesZeroMatrix() {
        Overlap o = new Overlap(3);
        assertEquals(3, o.getRowDimension());
        assertEquals(3, o.getColumnDimension());
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, o.getEntry(i, j), 1e-15);
            }
        }
    }
}
