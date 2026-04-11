package name.mjw.jquante.math.interpolater;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class ThreePointInterpolaterTest {

    private ThreePointInterpolater interpolater;

    @BeforeEach
    void setUp() {
        interpolater = new ThreePointInterpolater();
    }

    /**
     * x layout: [x_start, y_start, z_start, x_point, y_point, z_point, x_inc, y_inc, z_inc]
     * y layout: 8 corner values of the cube
     * delta[k] = (x_point[k] - x_start[k]) / x_inc[k]
     */

    @Test
    void constantFieldReturnsConstantValue() {
        double[] y = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
        double[] x = {0.0, 0.0, 0.0,  0.5, 0.5, 0.5,  1.0, 1.0, 1.0};
        assertEquals(5.0, interpolater.interpolate(y, x), 1e-10);
    }

    @Test
    void atLowerCornerReturnsFirstCornerValue() {
        // delta all zero: interpolation result = y[0] (lower corner value)
        double[] y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        double[] x = {0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  1.0, 1.0, 1.0};
        assertEquals(1.0, interpolater.interpolate(y, x), 1e-10);
    }

    @Test
    void linearGradientInXMidpoint() {
        // y[0..3]=0 (x=0 face), y[4..7]=2 (x=1 face); midpoint in x -> 1.0
        double[] y = {0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0};
        double[] x = {0.0, 0.0, 0.0,  0.5, 0.0, 0.0,  1.0, 1.0, 1.0};
        assertEquals(1.0, interpolater.interpolate(y, x), 1e-10);
    }

    @Test
    void linearGradientInXAtFullExtent() {
        // Point at x=1 face: should return y[4..7] (all = 2.0)
        double[] y = {0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0};
        double[] x = {0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  1.0, 1.0, 1.0};
        assertEquals(2.0, interpolater.interpolate(y, x), 1e-10);
    }

    @Test
    void numberOfRequiredYArgs() {
        assertEquals(8, interpolater.getNumberOfRequiredYArgs());
    }

    @Test
    void numberOfRequiredXArgs() {
        assertEquals(9, interpolater.getNumberOfRequiredXArgs());
    }
}
