package name.mjw.jquante.math.geom;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class BoundingBoxTest {

    private BoundingBox box;

    @BeforeEach
    void setUp() {
        box = new BoundingBox(new Vector3D(0, 0, 0), new Vector3D(2, 4, 6));
    }

    @Test
    void defaultConstructorIsZeroBox() {
        BoundingBox bb = new BoundingBox();
        assertEquals(0.0, bb.getXWidth(), 1e-10);
        assertEquals(0.0, bb.getYWidth(), 1e-10);
        assertEquals(0.0, bb.getZWidth(), 1e-10);
    }

    @Test
    void xWidth() {
        assertEquals(2.0, box.getXWidth(), 1e-10);
    }

    @Test
    void yWidth() {
        assertEquals(4.0, box.getYWidth(), 1e-10);
    }

    @Test
    void zWidth() {
        assertEquals(6.0, box.getZWidth(), 1e-10);
    }

    @Test
    void center() {
        Vector3D c = box.center();
        assertEquals(1.0, c.getX(), 1e-10);
        assertEquals(2.0, c.getY(), 1e-10);
        assertEquals(3.0, c.getZ(), 1e-10);
    }

    @Test
    void volume() {
        assertEquals(48.0, box.volume(), 1e-10);
    }

    @Test
    void totalSurfaceArea() {
        // 2*(2*4 + 4*6 + 6*2) = 2*(8 + 24 + 12) = 88
        assertEquals(88.0, box.totalSurfaceArea(), 1e-10);
    }

    @Test
    void containsInteriorPoint() {
        assertTrue(box.contains(new Vector3D(1, 2, 3)));
    }

    @Test
    void containsBoundaryPoint() {
        assertTrue(box.contains(new Vector3D(0, 0, 0)));
        assertTrue(box.contains(new Vector3D(2, 4, 6)));
    }

    @Test
    void doesNotContainExteriorPoint() {
        assertFalse(box.contains(new Vector3D(3, 2, 3)));
    }

    @Test
    void combineProducesEnclosingBox() {
        BoundingBox other = new BoundingBox(new Vector3D(-1, -1, -1), new Vector3D(1, 1, 1));
        BoundingBox combined = box.combine(other);
        assertEquals(-1.0, combined.getUpperLeft().getX(), 1e-10);
        assertEquals(6.0, combined.getBottomRight().getZ(), 1e-10);
    }

    @Test
    void expandSymmetric() {
        // expand(e) adds e on each side, so each width increases by 2*e
        BoundingBox expanded = box.expand(1.0);
        assertEquals(4.0, expanded.getXWidth(), 1e-10);
        assertEquals(6.0, expanded.getYWidth(), 1e-10);
        assertEquals(8.0, expanded.getZWidth(), 1e-10);
    }

    @Test
    void expandAsymmetric() {
        // expand(ex, ey, ez) adds ex on each x side, etc.
        BoundingBox expanded = box.expand(1.0, 2.0, 3.0);
        assertEquals(4.0, expanded.getXWidth(), 1e-10);
        assertEquals(8.0, expanded.getYWidth(), 1e-10);
        assertEquals(12.0, expanded.getZWidth(), 1e-10);
    }

    @Test
    void cornersHasEightPoints() {
        assertEquals(8, box.corners().length);
    }

    @Test
    void cornersFirstIsUpperLeft() {
        Vector3D[] corners = box.corners();
        assertEquals(0.0, corners[0].getX(), 1e-10);
        assertEquals(0.0, corners[0].getY(), 1e-10);
        assertEquals(0.0, corners[0].getZ(), 1e-10);
    }

    @Test
    void cloneProducesEqualBox() throws CloneNotSupportedException {
        BoundingBox cloned = (BoundingBox) box.clone();
        assertEquals(box.getXWidth(), cloned.getXWidth(), 1e-10);
        assertEquals(box.getYWidth(), cloned.getYWidth(), 1e-10);
        assertEquals(box.getZWidth(), cloned.getZWidth(), 1e-10);
    }

    @Test
    void toStringContainsCorners() {
        String s = box.toString();
        assertTrue(s.contains("UpperLeft"));
        assertTrue(s.contains("BottomRight"));
    }
}
