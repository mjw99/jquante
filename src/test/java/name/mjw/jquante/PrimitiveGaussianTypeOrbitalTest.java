package name.mjw.jquante;

import static org.junit.Assert.assertEquals;

import javax.vecmath.Point3d;

import org.junit.Test;

public class PrimitiveGaussianTypeOrbitalTest {

	private final double delta = 0.000001;

	@Test
	public void d() {
		PrimitiveGaussianTypeOrbital gto = new PrimitiveGaussianTypeOrbital(1.0);

		assertEquals(0.712705, gto.valueAtPoint(new Point3d(0, 0, 0)), delta);

	}

	@Test
	public void testTwo() {
		PrimitiveGaussianTypeOrbital gto = new PrimitiveGaussianTypeOrbital(
				1.0, new Point3d(0, 0, 0), 1, 0, 1);

		assertEquals(0.0, gto.valueAtPoint(new Point3d(0, 0, 0)), delta);

	}

}
