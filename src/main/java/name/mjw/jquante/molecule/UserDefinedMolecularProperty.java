package name.mjw.jquante.molecule;

import java.io.Serializable;

/**
 * A simple framework for adding arbitary user defined molecular property.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class UserDefinedMolecularProperty implements Serializable {

	/**
	 * Eclipse generated serialVersionUID
	 */
	private static final long serialVersionUID = -6598766569031535775L;

	/** Create a new instance of UserDefinedMolecularProperty */
	public UserDefinedMolecularProperty() {
	}

	/**
	 * Creates new instance of UserDefinedAtomProperty
	 * 
	 * @param name
	 *            Property name.
	 * @param value
	 *            Property value.
	 */
	public UserDefinedMolecularProperty(String name, Serializable value) {
		this.name = name;
		this.value = value;
	}

	/** The name identifying this user-defined molecular property. */
	protected String name;

	/**
	 * Get the value of name
	 * 
	 * @return the value of name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set the value of name
	 * 
	 * @param name
	 *            new value of name
	 */
	public void setName(String name) {
		this.name = name;
	}

	/** The value of this user-defined molecular property. */
	protected Serializable value;

	/**
	 * Get the value of value
	 * 
	 * @return the value of value
	 */
	public Serializable getValue() {
		return value;
	}

	/**
	 * Set the value of value
	 * 
	 * @param value
	 *            new value of value
	 */
	public void setValue(Serializable value) {
		this.value = value;
	}

	/** Whether this property depends on the machine or environment where the Molecule was created. */
	protected boolean originDependent;

	/**
	 * Get the value of originDependent
	 * 
	 * @return the value of originDependent
	 */
	public boolean isOriginDependent() {
		return originDependent;
	}

	/**
	 * Set the value of originDependent. (Note: origin dependent property is the
	 * one that is depended on the machine where Molecule object was created, an
	 * example is CommonUserDefinedMolecularProperty.SOURCE_FILE_NAME)
	 * 
	 * @param originDependent
	 *            new value of originDependent
	 */
	public void setOriginDependent(boolean originDependent) {
		this.originDependent = originDependent;
	}
}
