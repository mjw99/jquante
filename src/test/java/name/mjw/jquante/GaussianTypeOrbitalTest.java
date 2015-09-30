package name.mjw.jquante;

import static org.junit.Assert.*;

import javax.vecmath.Point3d;

import org.junit.Test;

public class GaussianTypeOrbitalTest {

	private final double delta = 0.000001;

	@Test
	public void testOne() {
		GaussianTypeOrbital gto = new GaussianTypeOrbital(1.0);

		assertEquals(0.712705, gto.valueAtPoint(new Point3d(0, 0, 0)), delta);

	}

	@Test
	public void testTwo() {
		GaussianTypeOrbital gto = new GaussianTypeOrbital(1.0, new Point3d(0,
				0, 0), 1, 0, 1);

		assertEquals(0.0, gto.valueAtPoint(new Point3d(0, 0, 0)), delta);

	}

}
