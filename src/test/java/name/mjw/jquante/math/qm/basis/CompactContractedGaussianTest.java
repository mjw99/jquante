package name.mjw.jquante.math.qm.basis;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

import org.junit.Test;

public class CompactContractedGaussianTest {

	private final double delta = 0.000001;

	@Test
	public void s() {
		ContractedGaussian cgto = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(0, 0, 0));

		cgto.addPrimitive(1.0, 1.0);
		cgto.normalize();

		CompactContractedGaussian ccgto = cgto.toCompactContractedGaussian();

		assertEquals(0, ccgto.l);
		assertEquals(0.0, ccgto.x, delta);
		assertEquals(1.0, ccgto.exponents[0], delta);

	}

	@Test
	public void sTO_3G_1S_H() {
		// https://bse.pnl.gov/bse/portal

		ContractedGaussian cgto = new ContractedGaussian(new Vector3D(0, 0, 0), new Power(0, 0, 0));

		cgto.addPrimitive(3.42525091, 0.154329);
		cgto.addPrimitive(0.62391373, 0.535328);
		cgto.addPrimitive(0.16885540, 0.444636);
		cgto.normalize();

		CompactContractedGaussian ccgto = cgto.toCompactContractedGaussian();

		assertEquals(3.42525091, ccgto.exponents[0], delta);
		assertEquals(0.62391373, ccgto.exponents[1], delta);
		assertEquals(0.16885540, ccgto.exponents[2], delta);
		
		assertEquals(0.154329, ccgto.coefficients[0], delta);
		assertEquals(0.535328, ccgto.coefficients[1], delta);
		assertEquals(0.444636, ccgto.coefficients[2], delta);


	}

}

