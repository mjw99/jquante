package name.mjw.jquante.molecule.event;

/**
 * The listeners who want to update them self if a concerned molecule object has
 * changed.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface MoleculeStateChangeListener extends java.util.EventListener {

	/**
	 * event fired when the state of the Molecule changes
	 * 
	 * @param event
	 *            - instance of MoleculeStateChangeEvent indicating the changes
	 *            in Molecule 's state.
	 */
	public void moleculeChanged(MoleculeStateChangeEvent event);

}