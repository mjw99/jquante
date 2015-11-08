package name.mjw.jquante.math;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MathUtilTest {

	@Test
	public void testFactorialZero() {

		assertEquals(1, MathUtil.factorial(0));
	}

	@Test
	public void testFactorialOne() {

		assertEquals(1, MathUtil.factorial(1));
	}

	@Test
	public void testFactorialTen() {

		assertEquals(3628800, MathUtil.factorial(10));
	}

	@Test
	public void testFactorial2MinusOne() {

		assertEquals(1, MathUtil.factorial2(-1));

	}

	@Test
	public void testFactorial2Zero() {

		assertEquals(1, MathUtil.factorial2(0));

	}

	@Test
	public void testFactorial2Nine() {

		assertEquals(945, MathUtil.factorial2(9));

	}

	@Test
	public void binomialFiveTwo() {
		assertEquals(10, MathUtil.binomial(5, 2));
	}

	@Test
	public void binomialTenFive() {
		assertEquals(252, MathUtil.binomial(10, 5));
	}

}
