package name.mjw.jquante.math.la;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.Vector;
import name.mjw.jquante.math.la.exception.SingularMatrixException;

public class GaussianEliminationTest {

	double delta = 0.000001;

	@Test
	public void testOne() {

		double[][] aValues = new double[][] { { 4 } };

		double[] bValues = new double[] { 8 };

		Matrix a = new Matrix(aValues);
		Vector b = new Vector(bValues);

		GaussianElimination gele = new GaussianElimination();
		gele.setMatrixA(a);
		gele.setVectorX(b);

		try {
			gele.findSolution();
			assertEquals(2.0, gele.getVectorX().getVector()[0], delta);

		} catch (SingularMatrixException e) {

			e.printStackTrace();
		}

	}

	@Test
	public void testTwo() {

		double[][] aValues = new double[][] { { 1, 1 }, { 2, 1 } };

		double[] bValues = new double[] { 10, 16 };

		Matrix a = new Matrix(aValues);
		Vector b = new Vector(bValues);

		GaussianElimination gele = new GaussianElimination();
		gele.setMatrixA(a);
		gele.setVectorX(b);

		try {
			gele.findSolution();
			assertEquals(6.0, gele.getVectorX().getVector()[0], delta);
			assertEquals(4.0, gele.getVectorX().getVector()[1], delta);

		} catch (SingularMatrixException e) {

			e.printStackTrace();
		}

	}

	@Test
	public void testThree() {

		double[][] aValues = new double[][] { { 1, 1, 1 }, { 2, 1, 2 }, { 1, 2, 3 } };

		double[] bValues = new double[] { 6, 10, 14 };

		Matrix a = new Matrix(aValues);
		Vector b = new Vector(bValues);

		GaussianElimination gele = new GaussianElimination();
		gele.setMatrixA(a);
		gele.setVectorX(b);

		try {
			gele.findSolution();
			assertEquals(1.0, gele.getVectorX().getVector()[0], delta);
			assertEquals(2.0, gele.getVectorX().getVector()[1], delta);
			assertEquals(2.0, gele.getVectorX().getVector()[1], delta);

		} catch (SingularMatrixException e) {

			e.printStackTrace();
		}

	}

	@Test
	public void testFour() {

		double[][] aValues = new double[][] { { 7, 1 }, { 5, 3 }, };

		double[] bValues = new double[] { 1, 1 };

		Matrix a = new Matrix(aValues);
		Vector b = new Vector(bValues);

		GaussianElimination gele = new GaussianElimination();
		gele.setMatrixA(a);
		gele.setVectorX(b);

		try {
			gele.findSolution();
			assertEquals(0.125, gele.getVectorX().getVector()[0], delta);
			assertEquals(0.125, gele.getVectorX().getVector()[1], delta);

		} catch (SingularMatrixException e) {

			e.printStackTrace();
		}

	}

}
