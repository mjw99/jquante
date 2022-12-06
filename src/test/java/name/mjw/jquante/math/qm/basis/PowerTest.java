package name.mjw.jquante.math.qm.basis;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

class PowerTest {

	private static Power s;
	private static Power p;
	private static Power d;
	private static Power f;

	@BeforeAll
	public static void setUp() {

		s = new Power(0, 0, 0);
		p = new Power(0, 1, 0);
		d = new Power(2, 0, 0);
		f = new Power(0, 2, 1);

	}

	@Test
	void testPowerS() {
		assertEquals(0, s.l());

	}

	@Test
	void testPowerP() {
		assertEquals(1, p.m());

	}

	@Test
	void testGetMaximumAngularMomentumD() {
		assertEquals(2, d.getMaximumAngularMomentum());

	}

	@Test
	void testGetTotalAngularMomentumF() {
		assertEquals(3, f.getTotalAngularMomentum());

	}

	@Test
	void testGetMinimumAngularMomentumF() {
		assertEquals(0, f.getMinimumAngularMomentum());

	}

	@Test
	void testToArrayF() {
		assertArrayEquals(new double[] { 0, 2, 1 }, f.toArray());

	}

}
