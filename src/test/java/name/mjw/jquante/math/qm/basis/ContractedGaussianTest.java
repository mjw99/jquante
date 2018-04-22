package name.mjw.jquante.math.qm.basis;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;

import org.junit.Test;

public class ContractedGaussianTest {

	private final double delta = 0.000001;

	@Test
	public void s() {
		ContractedGaussian cgto = new ContractedGaussian(new Vector3D(0, 0, 0),
				new Power(0, 0, 0));

		cgto.addPrimitive(1.0, 1.0);
		cgto.normalize();

		assertEquals(0.712705, cgto.amplitude(new Vector3D(0, 0, 0)), delta);

	}

	@Test
	public void sTO_3G_1S_H() {
		// https://bse.pnl.gov/bse/portal

		ContractedGaussian cgto = new ContractedGaussian(new Vector3D(0, 0, 0),
				new Power(0, 0, 0));

		cgto.addPrimitive(3.42525091, 0.154329);
		cgto.addPrimitive(0.62391373, 0.535328);
		cgto.addPrimitive(0.16885540, 0.444636);
		cgto.normalize();

		assertEquals(0.6282471373416881, cgto.amplitude(new Vector3D(0, 0, 0)),
				delta);

	}

}
