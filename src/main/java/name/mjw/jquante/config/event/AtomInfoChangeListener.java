package name.mjw.jquante.config.event;

import java.util.EventListener;

/**
 * The lisner interface, defining a mechanism for the observer to "listen" to
 * changes in AtomInfo object.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface AtomInfoChangeListener extends EventListener {

	/**
	 * The call back method for indicating the listeners of the said change.
	 * 
	 * @param event
	 *            - the event containig the required information.
	 */
	public void atomInfoChanged(AtomInfoChangeEvent event);

}