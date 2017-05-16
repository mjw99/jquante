package name.mjw.jquante.math.la.exception;

/**
 * Exception thrown when a matrix is singular, i.e. the matrix cannot be
 * inverted.
 * 
 * @author V. Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class SingularMatrixException extends Exception {

	/**
	 * Eclipse generated serialVersionUID
	 */
	private static final long serialVersionUID = 4254324891358597819L;

	/**
	 * Creates a new instance of <code>SingularMatrixException</code> without
	 * detail message.
	 */
	public SingularMatrixException() {
	}

	/**
	 * Constructs an instance of <code>SingularMatrixException</code> with the
	 * specified detail message.
	 * 
	 * @param msg
	 *            the detail message.
	 */
	public SingularMatrixException(String msg) {
		super(msg);
	}
}
