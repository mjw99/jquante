package name.mjw.jquante.math.qm.integral;

import static org.junit.Assert.assertEquals;

import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

import org.junit.Ignore;
import org.junit.Test;

public class RysTwoElectronTermTest {

	private final double delta = 0.000001;

	@Test
	@Ignore("To be implemented")
	public void testOne() {

		ContractedGaussian cgto1 = new ContractedGaussian(new Point3D(0, 0, 0),
				new Power(0, 0, 0));
		cgto1.addPrimitive(1.0, 1.0);
		cgto1.normalize();

		RysTwoElectronTerm e2 = new RysTwoElectronTerm();

		assertEquals(1.1283791633342477,
				e2.coulomb(cgto1, cgto1, cgto1, cgto1), delta);
	}

	@Test
	@Ignore("To be implemented")
	public void testTwo() {

		ContractedGaussian cgto1 = new ContractedGaussian(new Point3D(0, 0, 0),
				new Power(0, 0, 0));
		cgto1.addPrimitive(1.0, 1.0);
		cgto1.normalize();

		ContractedGaussian cgto2 = new ContractedGaussian(new Point3D(0, 0, 1),
				new Power(0, 0, 0));
		cgto2.addPrimitive(1.0, 1.0);
		cgto2.normalize();

		RysTwoElectronTerm e2 = new RysTwoElectronTerm();

		assertEquals(0.8427007900292194,
				e2.coulomb(cgto1, cgto1, cgto2, cgto2), delta);
	}
}
