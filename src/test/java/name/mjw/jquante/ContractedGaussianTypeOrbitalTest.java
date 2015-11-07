package name.mjw.jquante;

import static org.junit.Assert.*;

import javax.vecmath.Point3d;

import org.junit.Test;

public class ContractedGaussianTypeOrbitalTest {
	private final double delta = 0.000001;

	@Test
	public void testOne() {
		PrimitiveGaussianTypeOrbital pgto = new PrimitiveGaussianTypeOrbital(
				1.0);

		ContractedGaussianTypeOrbital cgto = new ContractedGaussianTypeOrbital();

		cgto.addPrimitiveGaussianTypeOrbital(1.0, pgto);

		assertEquals(0.712705, cgto.valueAtPoint(new Point3d(0, 0, 0)), delta);

	}

	@Test
	public void testSTO_3G_1S_H() {
		// https://bse.pnl.gov/bse/portal
		PrimitiveGaussianTypeOrbital pgto1 = new PrimitiveGaussianTypeOrbital(
				3.42525091);
		PrimitiveGaussianTypeOrbital pgto2 = new PrimitiveGaussianTypeOrbital(
				0.62391373);
		PrimitiveGaussianTypeOrbital pgto3 = new PrimitiveGaussianTypeOrbital(
				0.16885540);

		ContractedGaussianTypeOrbital cgto = new ContractedGaussianTypeOrbital();

		cgto.addPrimitiveGaussianTypeOrbital(0.154329, pgto1);
		cgto.addPrimitiveGaussianTypeOrbital(0.535328, pgto2);
		cgto.addPrimitiveGaussianTypeOrbital(0.444636, pgto3);

		// ~0.46?
		// assertEquals(0.712705, cgto.valueAtPoint(new Point3d(0, 0, 0)),
		// delta);

	}
}
