package name.mjw.jquante;

public class Utils {

	public static int factorial(int x) {

		int fact = 1;
		while (x > 1) {
			fact = fact * x;
			x = x - 1;
		}

		return fact;
	}

	public static int doubleFactorial(int x) {

		int fact = 1;
		while (x > 1) {
			fact = fact * x;
			x = x - 2;
		}
		return fact;
	}

	/**
	 * Binomial coefficient
	 * 
	 * Number of ways to choose k elements from a set of n elements.
	 * 
	 * @param n
	 *            number of n elements in set
	 * @param k
	 *            number of k elements in set
	 * @return binomial coefficient
	 * 
	 *         https://en.wikipedia.org/wiki/Binomial_coefficient
	 */
	public static int binomial(int n, int k) {
		if (n == k) {
			return 1;
		}

		return factorial(n) / (factorial(k) * factorial(n - k));
	}

}
