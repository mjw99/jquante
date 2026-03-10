package name.mjw.jquante.math.interpolater;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class CubicInterpolaterTest {

    private CubicInterpolater interpolater;

    @BeforeEach
    void setUp() {
        interpolater = new CubicInterpolater();
    }

    @Test
    void interpolateAtZeroReturnsY1() {
        // mu=0: mu2=0, result = a1*0 + a2*0 + y[1] = y[1]
        double[] y = {1.0, 2.0, 3.0, 4.0};
        assertEquals(2.0, interpolater.interpolate(y, new double[]{0.0}), 1e-10);
    }

    @Test
    void interpolateAtOneMidpointForLinearData() {
        // For purely linear data y={0,1,2,3}: a0=3-2-0+1=2, a1=0-1-2=-3, a2=2-0=2, a3=1
        // at mu=0.5: mu2=0.25
        // result = 2*0.5*0.25 + (-3)*0.25 + 2*0.5 + 1 = 0.25 - 0.75 + 1 + 1 = 1.5
        double[] y = {0.0, 1.0, 2.0, 3.0};
        assertEquals(1.5, interpolater.interpolate(y, new double[]{0.5}), 1e-10);
    }

    @Test
    void interpolateConstantDataReturnsConstant() {
        // For y={5,5,5,5}: a0=5-5-5+5=0, a1=5-5-0=0, a2=5-5=0, a3=5
        // result = 5 for any mu
        double[] y = {5.0, 5.0, 5.0, 5.0};
        assertEquals(5.0, interpolater.interpolate(y, new double[]{0.3}), 1e-10);
        assertEquals(5.0, interpolater.interpolate(y, new double[]{0.7}), 1e-10);
    }

    @Test
    void numberOfRequiredYArgs() {
        assertEquals(4, interpolater.getNumberOfRequiredYArgs());
    }

    @Test
    void numberOfRequiredXArgs() {
        assertEquals(1, interpolater.getNumberOfRequiredXArgs());
    }

    @Test
    void toStringDescribesType() {
        assertEquals("Cubic interpolater", interpolater.toString());
    }
}
