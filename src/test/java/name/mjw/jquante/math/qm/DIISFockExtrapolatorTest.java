package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class DIISFockExtrapolatorTest {

    private DIISFockExtrapolator diis;
    // Identity overlap so that FPS - SPF = F*P - P*F; with zero density this is 0
    private Overlap identityOverlap;
    // Zero density: error vector = F*0*S - S*0*F = 0, so isDiisStarted stays false (bootstrap mode)
    private Density zeroDensity;

    @BeforeEach
    void setUp() {
        diis = new DIISFockExtrapolator();
        identityOverlap = new Overlap(2);
        identityOverlap.setEntry(0, 0, 1.0);
        identityOverlap.setEntry(1, 1, 1.0);
        zeroDensity = new Density(2);
    }

    @Test
    void defaultErrorThreshold() {
        assertEquals(0.00001, diis.getErrorThreshold(), 1e-10);
    }

    @Test
    void setAndGetErrorThreshold() {
        diis.setErrorThreshold(0.001);
        assertEquals(0.001, diis.getErrorThreshold(), 1e-10);
    }

    @Test
    void bootstrapFirstCallReturnsCurrentFock() {
        Fock f1 = new Fock(new double[][]{{2.0, 0.5}, {0.5, 3.0}});
        Fock result = diis.next(f1, identityOverlap, zeroDensity);

        // First call with oldFock=null: current Fock is returned unchanged
        assertEquals(2.0, result.getEntry(0, 0), 1e-10);
        assertEquals(0.5, result.getEntry(0, 1), 1e-10);
        assertEquals(3.0, result.getEntry(1, 1), 1e-10);
    }

    @Test
    void bootstrapSecondCallReturnsAverageOfOldAndCurrent() {
        Fock f1 = new Fock(new double[][]{{2.0, 0.0}, {0.0, 2.0}});
        Fock f2 = new Fock(new double[][]{{4.0, 0.0}, {0.0, 4.0}});

        diis.next(f1, identityOverlap, zeroDensity);   // sets oldFock = f1
        Fock result = diis.next(f2, identityOverlap, zeroDensity);  // returns 0.5*f1 + 0.5*f2

        assertEquals(3.0, result.getEntry(0, 0), 1e-10);
        assertEquals(3.0, result.getEntry(1, 1), 1e-10);
        assertEquals(0.0, result.getEntry(0, 1), 1e-10);
    }

    @Test
    void bootstrapSecondCallWithAsymmetricFocks() {
        Fock f1 = new Fock(new double[][]{{0.0, 0.0}, {0.0, 0.0}});
        Fock f2 = new Fock(new double[][]{{6.0, 2.0}, {2.0, 4.0}});

        diis.next(f1, identityOverlap, zeroDensity);
        Fock result = diis.next(f2, identityOverlap, zeroDensity);

        // 0.5*[[0,0],[0,0]] + 0.5*[[6,2],[2,4]] = [[3,1],[1,2]]
        assertEquals(3.0, result.getEntry(0, 0), 1e-10);
        assertEquals(1.0, result.getEntry(0, 1), 1e-10);
        assertEquals(1.0, result.getEntry(1, 0), 1e-10);
        assertEquals(2.0, result.getEntry(1, 1), 1e-10);
    }

    @Test
    void diisStepStartsAtZero() {
        assertEquals(0, diis.diisStep);
    }

    @Test
    void bootstrapDoesNotIncrementDiisStep() {
        Fock f = new Fock(new double[][]{{1.0, 0.0}, {0.0, 1.0}});
        diis.next(f, identityOverlap, zeroDensity);
        diis.next(f, identityOverlap, zeroDensity);
        // Bootstrap mode never increments diisStep
        assertEquals(0, diis.diisStep);
    }
}
