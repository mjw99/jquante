package name.mjw.jquante.math.qm.integral;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import name.mjw.jquante.math.geom.Point3D;

public class IntegralsUtilTest {

	private final double delta = 0.00001;

	@Test
	public void testGaussianProductCenterOne() {

		assertEquals(new Point3D(0.0, 0.0, 0.0),
				IntegralsUtil.gaussianProductCenter(1.0, new Point3D(0, 0, 0), 1.0, new Point3D(0, 0, 0)));
	}

	@Test
	public void testGaussianProductCenterTwo() {

		assertEquals(new Point3D(0.5, 0.5, 0.5),
				IntegralsUtil.gaussianProductCenter(1.0, new Point3D(0, 0, 0), 1.0, new Point3D(1, 1, 1)));
	}

	@Test
	public void testIjkl2intindex0000() {

		assertEquals(0, IntegralsUtil.ijkl2intindex(0, 0, 0, 0));
	}

	@Test
	public void testIjkl2intindex1000() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(1, 0, 0, 0));
	}

	@Test
	public void testIjkl2intindex0100() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(0, 1, 0, 0));
	}

	@Test
	public void testIjkl2intindex0010() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(0, 0, 1, 0));
	}

	@Test
	public void testIjkl2intindex0001() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(0, 0, 0, 1));
	}

	@Test
	public void testIjkl2intindex1001() {

		assertEquals(2, IntegralsUtil.ijkl2intindex(1, 0, 0, 1));
	}

	@Test
	public void testIjkl2intindex1100() {

		assertEquals(3, IntegralsUtil.ijkl2intindex(1, 1, 0, 0));
	}

	@Test
	public void testIjkl2intindex0011() {

		assertEquals(3, IntegralsUtil.ijkl2intindex(0, 0, 1, 1));
	}

	@Test
	public void computeFGamma() {
		assertEquals(1.0, IntegralsUtil.computeFGamma(0, 0), delta);

	}

	@Test
	public void gammaIncompleteOne() {
		assertEquals(1.49365, IntegralsUtil.gammaIncomplete(0.5, 1), delta);

	}

	@Test
	public void gammaIncompleteTwo() {
		assertEquals(0.6545103, IntegralsUtil.gammaIncomplete(1.5, 2), delta);

	}

	@Test
	public void gammaIncompleteThree() {
		assertEquals(0, IntegralsUtil.gammaIncomplete(2.5, 1E-12), delta);

	}

}
