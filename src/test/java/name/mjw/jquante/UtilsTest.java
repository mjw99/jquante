package name.mjw.jquante;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class UtilsTest {

	@Test
	public void testFactorialZero() {

		assertEquals(1, Utils.factorial(1));
	}

	@Test
	public void testFactorialOne() {

		assertEquals(1, Utils.factorial(1));
	}

	@Test
	public void testFactorialTen() {

		assertEquals(3628800, Utils.factorial(10));
	}

	@Test
	public void testDoubleFactorialMinusOne() {

		assertEquals(1, Utils.doubleFactorial(-1));

	}

	@Test
	public void testDoubleFactorialZero() {

		assertEquals(1, Utils.doubleFactorial(0));

	}

	@Test
	public void testDoubleFactorialNine() {

		assertEquals(945, Utils.doubleFactorial(9));

	}

	@Test
	public void binomialFiveTwo() {
		assertEquals(10, Utils.binomial(5, 2));
	}

	@Test
	public void binomialTenFive() {
		assertEquals(252, Utils.binomial(10, 5));
	}

}
