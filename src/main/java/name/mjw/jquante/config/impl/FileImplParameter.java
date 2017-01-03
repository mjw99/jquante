package name.mjw.jquante.config.impl;

import java.io.File;

import name.mjw.jquante.config.Parameter;

/**
 * A simple implementation of Parameter interface for storing most rescently
 * opened workspace files.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class FileImplParameter implements Parameter {

	private File theFile;

	/**
	 * Creates a new instance of FileImplParameter
	 * 
	 * @param theFile
	 *            Name of the File.
	 */
	public FileImplParameter(File theFile) {
		this.theFile = theFile;
	}

	/**
	 * Getter for property value.
	 * 
	 * @return Value of property value.
	 */
	@Override
	public Object getValue() {
		return theFile;
	}

	/**
	 * Setter for property value.
	 * 
	 * @param value
	 *            New value of property value.
	 * 
	 */
	@Override
	public void setValue(Object value) {
		theFile = (File) value;
	}

}