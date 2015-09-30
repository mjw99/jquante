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

	static public int doubleFactorial(int x) {

		int fact = 1;
		while (x > 1) {
			fact = fact * x;
			x = x - 2;
		}
		return fact;
	}

}
