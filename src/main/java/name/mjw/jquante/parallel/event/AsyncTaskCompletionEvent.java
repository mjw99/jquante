package name.mjw.jquante.parallel.event;

/**
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class AsyncTaskCompletionEvent extends java.util.EventObject {

	/**
	 * Eclipse generated serialVersionUID
	 */
	private static final long serialVersionUID = 1028970494565784834L;

	/**
	 * Creates a new instance of AsyncTaskCompletionEvent
	 * 
	 * @param source
	 *            The object upon which the Event in question initially occurred
	 *            upon.
	 */
	public AsyncTaskCompletionEvent(Object source) {
		super(source);
	}

	protected Object returnValue;

	/**
	 * Get the value of returnValue
	 * 
	 * @return the value of returnValue
	 */
	public Object getReturnValue() {
		return returnValue;
	}

	/**
	 * Set the value of returnValue
	 * 
	 * @param returnValue
	 *            new value of returnValue
	 */
	public void setReturnValue(Object returnValue) {
		this.returnValue = returnValue;
	}
}
