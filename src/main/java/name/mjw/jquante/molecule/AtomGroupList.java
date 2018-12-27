package name.mjw.jquante.molecule;

import java.util.ArrayList;
import java.util.Iterator;


/**
 * An interface defining a list of AtomGroup instances.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface AtomGroupList {

	/**
	 * Add a new atom group.
	 * 
	 * @param atomGroup
	 *            to be added to the group list
	 */
	public void addGroup(AtomGroup atomGroup);

	/**
	 * Removes an atom group from the current list.
	 * 
	 * @param atomGroup
	 *            to be removes from the group list
	 */
	public void removeGroup(AtomGroup atomGroup);

	/**
	 * get the list of all group as an iterator object
	 * 
	 * @return Iterator object containing instance of AtomGroup
	 */
	public Iterator<AtomGroup> getGroups();

	/**
	 * get the size of this list
	 * 
	 * @return the size of this list
	 */
	public int getSize();

	/**
	 * Getter for property atomList.
	 * 
	 * @return Value of property atomList.
	 */
	public ArrayList<AtomGroup> getGroupList();

}