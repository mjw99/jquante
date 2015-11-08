package name.mjw.jquante.math.qm.basis;

import static org.junit.Assert.*;

import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PrimitiveGaussian;

import org.junit.Test;

public class PrimitiveGaussianTest {

	private final double delta = 0.000001;

	@Test
	public void testOne() {
		PrimitiveGaussian gto = new PrimitiveGaussian(new Point3D(0, 0, 0),
				new Power(0, 0, 0), 1.0, 1.0);

		assertEquals(0.712705, gto.amplitude(new Point3D(0, 0, 0)), delta);

	}

	@Test
	public void testTwo() {
		PrimitiveGaussian gto = new PrimitiveGaussian(new Point3D(0, 0, 0),
				new Power(1, 0, 1), 1.0, 1.0);

		assertEquals(0.0, gto.amplitude(new Point3D(0, 0, 0)), delta);

	}
}
