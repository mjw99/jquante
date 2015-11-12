package name.mjw.jquante.math.qm.integral;

import static org.junit.Assert.*;

import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

import org.junit.Test;

public class HGPTwoElectronTermTest {
	private final double delta = 0.000001;

	@Test
	public void testOne() {

		ContractedGaussian cgto1 = new ContractedGaussian(new Point3D(0, 0, 0),
				new Power(0, 0, 0));
		cgto1.addPrimitive(1.0, 1.0);
		cgto1.normalize();

		HGPTwoElectronTerm e2 = new HGPTwoElectronTerm();

		assertEquals(1.1283791633342477,
				e2.coulomb(cgto1, cgto1, cgto1, cgto1), delta);
	}

	@Test
	public void testTwo() {

		ContractedGaussian cgto1 = new ContractedGaussian(new Point3D(0, 0, 0),
				new Power(0, 0, 0));
		cgto1.addPrimitive(1.0, 1.0);
		cgto1.normalize();

		ContractedGaussian cgto2 = new ContractedGaussian(new Point3D(0, 0, 1),
				new Power(0, 0, 0));
		cgto2.addPrimitive(1.0, 1.0);
		cgto2.normalize();

		HGPTwoElectronTerm e2 = new HGPTwoElectronTerm();

		assertEquals(0.8427007900292194,
				e2.coulomb(cgto1, cgto1, cgto2, cgto2), delta);
	}
}
