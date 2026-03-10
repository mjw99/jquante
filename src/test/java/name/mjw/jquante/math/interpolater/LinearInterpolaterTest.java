package name.mjw.jquante.math.interpolater;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class LinearInterpolaterTest {

    private LinearInterpolater interpolater;

    @BeforeEach
    void setUp() {
        interpolater = new LinearInterpolater();
    }

    @Test
    void interpolateAtZeroReturnFirstValue() {
        // mu=0 => y[0]*(1-0) + y[1]*0 = y[0]
        assertEquals(3.0, interpolater.interpolate(new double[]{3.0, 7.0}, new double[]{0.0}), 1e-10);
    }

    @Test
    void interpolateAtOneReturnsSecondValue() {
        // mu=1 => y[0]*0 + y[1]*1 = y[1]
        assertEquals(7.0, interpolater.interpolate(new double[]{3.0, 7.0}, new double[]{1.0}), 1e-10);
    }

    @Test
    void interpolateAtMidpoint() {
        // mu=0.5 => (3+7)/2 = 5
        assertEquals(5.0, interpolater.interpolate(new double[]{3.0, 7.0}, new double[]{0.5}), 1e-10);
    }

    @Test
    void interpolateAtQuarter() {
        // mu=0.25 => 3*0.75 + 7*0.25 = 2.25 + 1.75 = 4.0
        assertEquals(4.0, interpolater.interpolate(new double[]{3.0, 7.0}, new double[]{0.25}), 1e-10);
    }

    @Test
    void numberOfRequiredYArgs() {
        assertEquals(2, interpolater.getNumberOfRequiredYArgs());
    }

    @Test
    void numberOfRequiredXArgs() {
        assertEquals(1, interpolater.getNumberOfRequiredXArgs());
    }

    @Test
    void toStringDescribesType() {
        assertEquals("Linear interpolater", interpolater.toString());
    }
}
