package name.mjw.jquante.math.interpolater;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class TriLinearInterpolaterTest {

    private TriLinearInterpolater interpolater;

    @BeforeEach
    void setUp() {
        interpolater = new TriLinearInterpolater();
    }

    /**
     * x array layout: x1,y1,z1, x0,y0,z0, incX,incY,incZ
     * y array: Vx1y1z1, Vx1y1z2, Vx1y2z1, Vx1y2z2, Vx2y1z1, Vx2y1z2, Vx2y2z1, Vx2y2z2
     */
    @Test
    void constantGridReturnsConstantValue() {
        // All corner values = 5, point exactly at lower corner, increments = 1
        double[] y = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
        // x1=0,y1=0,z1=0, x0=0,y0=0,z0=0, incX=1,incY=1,incZ=1 => deltas all = 0
        double[] x = {0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  1.0, 1.0, 1.0};
        assertEquals(5.0, interpolater.interpolate(y, x), 1e-10);
    }

    @Test
    void midpointOfLinearGradientInX() {
        // Values vary linearly in x: lower face = 0, upper face = 2
        // y: [Vx1y1z1=0, Vx1y1z2=0, Vx1y2z1=0, Vx1y2z2=0, Vx2y1z1=2, Vx2y1z2=2, Vx2y2z1=2, Vx2y2z2=2]
        double[] y = {0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0};
        // x1=0, y1=0, z1=0, x0=0.5, y0=0, z0=0, incX=1, incY=1, incZ=1
        // delta[0] = (0.5-0)/1 = 0.5
        double[] x = {0.0, 0.0, 0.0,  0.5, 0.0, 0.0,  1.0, 1.0, 1.0};
        assertEquals(1.0, interpolater.interpolate(y, x), 1e-10);
    }

    @Test
    void numberOfRequiredYArgs() {
        assertEquals(8, interpolater.getNumberOfRequiredYArgs());
    }

    @Test
    void numberOfRequiredXArgs() {
        assertEquals(9, interpolater.getNumberOfRequiredXArgs());
    }

    @Test
    void toStringDescribesType() {
        assertEquals("Tri linear interpolater", interpolater.toString());
    }
}
