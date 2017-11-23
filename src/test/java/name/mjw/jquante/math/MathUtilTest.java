package name.mjw.jquante.math;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MathUtilTest {
	double diff = 0.00001;

	@Test
	public void factorialZero() {

		assertEquals(1, MathUtil.factorial(0));
	}

	@Test
	public void factorialOne() {

		assertEquals(1, MathUtil.factorial(1));
	}

	@Test
	public void factorialTen() {

		assertEquals(3628800, MathUtil.factorial(10));
	}

	@Test
	public void factorial2MinusOne() {

		assertEquals(1, MathUtil.factorial2(-1));

	}

	@Test
	public void factorial2Zero() {

		assertEquals(1, MathUtil.factorial2(0));

	}

	@Test
	public void factorial2Nine() {

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

	@Test
	public void degreesToRadians() {
		assertEquals(90.0, MathUtil.toDegrees(Math.PI / 2), diff);

	}

	@Test
	public void radiansToDegrees() {
		assertEquals(Math.PI / 4, MathUtil.toRadians(45), diff);

	}

}
