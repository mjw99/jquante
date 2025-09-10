package name.mjw.jquante.math.qm.basis;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

public class OrbitalTest {

	private final double delta = 0.000001;
	private static Orbital orbital;

	private static Orbital largerOrbital;

	@BeforeAll
	public static void setUp() {
		orbital = new Orbital("S");
		orbital.addCoefficient(1.0);
		orbital.addExponent(2.0);

		largerOrbital = new Orbital("P");
		largerOrbital.addEntry(1.0, 2.0);
		largerOrbital.addEntry(3.0, 4.0);

	}

	@Test
	void testGetType() {
		assertEquals("S", orbital.getType());

	}

	@Test
	void testGetCoefficients() {
		assertEquals(1, orbital.getCoefficients().size());
		assertEquals(1.0, orbital.getCoefficients().get(0), delta);

	}

	@Test
	void testGetExponent() {
		assertEquals(1, orbital.getExponents().size());
		assertEquals(2.0, orbital.getExponents().get(0), delta);

	}

	@Test
	void testAddEntry() {
		assertEquals(2, largerOrbital.getCoefficients().size());

		assertEquals(1.0, largerOrbital.getCoefficients().get(0), delta);
		assertEquals(2.0, largerOrbital.getExponents().get(0), delta);

		assertEquals(3.0, largerOrbital.getCoefficients().get(1), delta);
		assertEquals(4.0, largerOrbital.getExponents().get(1), delta);

	}

}
