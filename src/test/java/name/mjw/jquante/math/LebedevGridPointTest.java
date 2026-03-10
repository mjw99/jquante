package name.mjw.jquante.math;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class LebedevGridPointTest {

    private static final double DELTA = 1e-10;

    // A point on the unit sphere: (1, 0, 0) with weight 0.5
    private LebedevGridPoint point;

    @BeforeEach
    void setUp() {
        point = new LebedevGridPoint(1.0, 0.0, 0.0, 0.5);
    }

    @Test
    void getX() {
        assertEquals(1.0, point.x(), DELTA);
    }

    @Test
    void getY() {
        assertEquals(0.0, point.y(), DELTA);
    }

    @Test
    void getZ() {
        assertEquals(0.0, point.z(), DELTA);
    }

    @Test
    void getWeight() {
        assertEquals(0.5, point.weight(), DELTA);
    }

    @Test
    void getVector3D() {
        Vector3D v = point.getVector3D();
        assertEquals(1.0, v.getX(), DELTA);
        assertEquals(0.0, v.getY(), DELTA);
        assertEquals(0.0, v.getZ(), DELTA);
    }

    @Test
    void getXYZ() {
        RealVector v = point.getXYZ();
        assertEquals(3, v.getDimension());
        assertEquals(1.0, v.getEntry(0), DELTA);
        assertEquals(0.0, v.getEntry(1), DELTA);
        assertEquals(0.0, v.getEntry(2), DELTA);
    }

    @Test
    void getPhiForNorthPole() {
        // (0, 0, 1): phi = acos(1) = 0
        LebedevGridPoint north = new LebedevGridPoint(0.0, 0.0, 1.0, 1.0);
        assertEquals(0.0, north.getPhi(), DELTA);
    }

    @Test
    void getPhiForEquator() {
        // z=0: phi = acos(0) = PI/2
        LebedevGridPoint equator = new LebedevGridPoint(1.0, 0.0, 0.0, 1.0);
        assertEquals(Math.PI / 2.0, equator.getPhi(), DELTA);
    }

    @Test
    void getThetaForXAxis() {
        // (1, 0, 0): theta = atan2(0, 1) = 0
        assertEquals(0.0, point.getTheta(), DELTA);
    }

    @Test
    void getThetaForYAxis() {
        // (0, 1, 0): theta = atan2(1, 0) = PI/2
        LebedevGridPoint yAxis = new LebedevGridPoint(0.0, 1.0, 0.0, 1.0);
        assertEquals(Math.PI / 2.0, yAxis.getTheta(), DELTA);
    }
}
