package name.mjw.jquante.math.qm.basis;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

public class OrbitalTest {

	private final double delta = 0.000001;
	private static Orbital orbital;

	@BeforeAll
	public static void setUp() {
		orbital = new Orbital("S");
		orbital.addCoefficient(1.0);
		orbital.addExponent(2.0);
	}

	@Test
	void testGetType() {
		assertEquals("S", orbital.getType());

	}

	@Test
	void testGetCoefficients() {
		assertEquals(1.0, orbital.getCoefficients().get(0), delta);

	}

	@Test
	void testGetExponent() {
		assertEquals(2.0, orbital.getExponents().get(0), delta);

	}

}
