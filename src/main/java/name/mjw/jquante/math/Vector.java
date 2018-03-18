package name.mjw.jquante.math;

import java.util.stream.IntStream;

/**
 * Defines an R^n space vector.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class Vector implements Cloneable {

	/** Holds value of property vector. */
	protected double[] data;

	/**
	 * Creates a new instance of Vector
	 * 
	 * @param size
	 *            the vector.length of this vector
	 */
	public Vector(int size) {
		data = new double[size];
	}

	/**
	 * Creates a new instance of Vector from an double array
	 * 
	 * @param array
	 *            the array from which this vector will be made
	 */
	public Vector(double[] array) {
		this.data = array;
	}

	/**
	 * Creates a new instance of Vector from a square Matrix
	 * 
	 * @param a
	 *            the square matrix from which this vector will be made
	 */
	public Vector(Matrix a) {
		double[][] matrix = a.getData();

		this.data = new double[matrix.length * matrix.length];

		int ii = 0;
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				this.data[ii] = matrix[i][j];
				ii++;
			}
		}
	}

	/**
	 * The size of this vector
	 * 
	 * @return the size of this vector
	 * 
	 */
	public int getSize() {
		return this.data.length;
	}

	/**
	 * Indexed getter for property vector.
	 * 
	 * @param index
	 *            Index of the property.
	 * @return Value of the property at <CODE>index</CODE>.
	 * 
	 */
	public double getVector(int index) {
		return this.data[index];
	}

	/**
	 * Getter for property vector.
	 * 
	 * @return Value of property vector.
	 * 
	 */
	public double[] getVector() {
		return this.data;
	}

	/**
	 * Indexed setter for property vector.
	 * 
	 * @param index
	 *            Index of the property.
	 * @param vector
	 *            New value of the property at <CODE>index</CODE>.
	 * 
	 */
	public void setVector(int index, double vector) {
		this.data[index] = vector;
	}

	/**
	 * Setter for property vector.
	 * 
	 * @param vector
	 *            New value of property vector.
	 * 
	 */
	public void setVector(double[] vector) {
		this.data = vector;
	}

	/**
	 * Make this vector a Null vector
	 */
	public void makeZero() {
		IntStream.range(0, data.length).parallel().forEach(id -> data[id] = 0.0);
	}

	/**
	 * addition of two vectors
	 * 
	 * @param b
	 *            Vector to be added to the current vector
	 * @return the addition of two vector or null if that is not possible
	 */
	public Vector add(Vector b) {
		if (this.data.length != b.data.length) {
			return null;
		}

		Vector result = new Vector(data.length);

		IntStream.range(0, result.data.length).parallel()
				.forEach(id -> result.data[id] = this.data[id] + b.data[id]);

		return result;
	}

	/**
	 * Subtraction of two vectors (this - b)
	 * 
	 * @param b
	 *            the vector to be subtracted from the current vector
	 * @return the addition of two vector or null if that is not possible
	 */
	public Vector sub(Vector b) {
		if (this.data.length != b.data.length) {
			return null;
		}

		Vector result = new Vector(data.length);

		IntStream.range(0, result.data.length).parallel()
				.forEach(id -> result.data[id] = this.data[id] - b.data[id]);

		return result;
	}

	/**
	 * multiplication of this vector with a scalar k
	 * 
	 * @param k
	 *            the scalar to be multiplied to this vector
	 * @return the result !
	 */
	public Vector mul(double k) {
		Vector result = new Vector(data.length);

		IntStream.range(0, result.data.length).parallel().forEach(id -> result.data[id] = this.data[id] * k);

		return result;
	}

	/**
	 * The dot product of two vectors (a.b)
	 * 
	 * @param b
	 *            the vector with which the dot product is to be evaluated
	 * @return a double value which is the result of the dot product
	 */
	public double dot(Vector b) {
		if (this.data.length != b.data.length) {
			return Double.NaN;
		}
		return IntStream.range(0, data.length).parallel().mapToDouble(id -> this.data[id] * b.data[id]).reduce(0,
				Double::sum);
	}

	/**
	 * The magnitude of this vector
	 * 
	 * @return the magnitude (length) of this vector.
	 */
	public double magnitude() {
		double length;

		length = IntStream.range(0, data.length).parallel().mapToDouble(id -> this.data[id] * this.data[id])
				.reduce(0, Double::sum);

		return Math.sqrt(length);
	}

	/**
	 * get the normalised form of this vector
	 * 
	 * @return the normalised form of this vector
	 */
	public Vector normalize() {
		double magnitude = magnitude();

		Vector n = new Vector(data.length);

		IntStream.range(0, data.length).parallel().mapToDouble(id -> n.data[id] = this.data[id] / magnitude);

		return n;
	}

	/**
	 * negate this vector
	 * 
	 * @return the reverse vector
	 */
	public Vector negate() {
		Vector n = new Vector(data.length);

		IntStream.range(0, data.length).parallel().mapToDouble(id -> n.data[id] = -this.data[id]);

		return n;
	}

	/**
	 * clone this vector ;) Cloning is getting interesting! :)
	 * 
	 * @return the clone
	 */
	@Override
	public Object clone() {
		Vector theCopy = new Vector(data.length);

		for (int i = 0; i < data.length; i++) {
			theCopy.data[i] = data[i];
		}

		return theCopy;
	}

	/**
	 * the overridden toString() method
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < data.length; i++) {
			sb.append(data[i]);
			sb.append("\n");
		}

		return sb.toString();
	}

	/**
	 * Returns the max norm of this matrix
	 * 
	 * @return max(abs(matrix))
	 */
	public double maxNorm() {

		return IntStream.range(0, data.length).parallel().mapToDouble(id -> this.data[id]).reduce(0, Double::max);

	}
}