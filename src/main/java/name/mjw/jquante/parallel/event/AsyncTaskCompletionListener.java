package name.mjw.jquante.parallel.event;

/**
 * Listener mechanism for handling async task completion.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface AsyncTaskCompletionListener extends java.util.EventListener {

	/**
	 * Called when an registered async task completes
	 * 
	 * @param atce
	 *            the AsyncTaskCompletionEvent object
	 */
	public void asyncTaskCompleted(AsyncTaskCompletionEvent atce);
}
