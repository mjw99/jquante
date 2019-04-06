package name.mjw.jquante.config.impl;

import java.awt.Color;
import java.beans.PropertyVetoException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashMap;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import name.mjw.jquante.common.EventListenerList;
import name.mjw.jquante.common.Utility;
import name.mjw.jquante.common.resource.StringResource;
import name.mjw.jquante.config.Configuration;
import name.mjw.jquante.config.Parameter;
import name.mjw.jquante.config.event.AtomInfoChangeEvent;
import name.mjw.jquante.config.event.AtomInfoChangeListener;

/**
 * The default AtomProperty configuration.
 * 
 * .. follows a singleton pattern. .. and an observer pattern for notifying the
 * registered classes of the changes at runtime.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class AtomInfo implements Configuration {

	private static AtomInfo atomInfo;

	private static final int DEFAULT_TABLE_SIZE = 90;

	private String element = "element";
	private String symbol = "symbol";

	/** Holds value of property nameTable. */
	private HashMap<String, String> nameTable;
	private HashMap<String, String> originalNameTable;

	/** Holds value of property atomicNumberTable. */
	private HashMap<String, Integer> atomicNumberTable;
	private HashMap<String, Integer> originalAtomicNumberTable;

	/** Holds value of property atomicWeightTable. */
	private HashMap<String, Double> atomicWeightTable;
	private HashMap<String, Double> originalAtomicWeightTable;

	/** Holds value of property covalentRadiusTable. */
	private HashMap<String, Double> covalentRadiusTable;
	private HashMap<String, Double> originalCovalentRadiusTable;

	/** Holds value of property vdwRadiusTable. */
	private HashMap<String, Double> vdwRadiusTable;
	private HashMap<String, Double> originalVdwRadiusTable;

	/** Holds value of property defaultValencyTable. */
	private HashMap<String, Integer> defaultValencyTable;
	private HashMap<String, Integer> originalDefaultValencyTable;

	/** Holds value of property weakBondAngleTable. */
	private HashMap<String, Double> weakBondAngleTable;
	private HashMap<String, Double> originalWeakBondAngleTable;

	/** Holds value of property colorTable. */
	private HashMap<String, Color> colorTable;
	private HashMap<String, Color> originalColorTable;

	/** Holds value of property doubleBondOverlapTable. */
	private HashMap<String, Double> doubleBondOverlapTable;
	private HashMap<String, Double> originalDoubleBondOverlapTable;

	/** Utility field used by event firing mechanism. */
	private EventListenerList<AtomInfoChangeListener> listenerList = null;

	/** Creates a new instance of AtomInfo */
	public AtomInfo() {
		nameTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		atomicNumberTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		atomicWeightTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		covalentRadiusTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		vdwRadiusTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		weakBondAngleTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		colorTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		defaultValencyTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		doubleBondOverlapTable = new HashMap<>(DEFAULT_TABLE_SIZE);

		originalNameTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalAtomicNumberTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalAtomicWeightTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalCovalentRadiusTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalVdwRadiusTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalWeakBondAngleTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalColorTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalDefaultValencyTable = new HashMap<>(DEFAULT_TABLE_SIZE);
		originalDoubleBondOverlapTable = new HashMap<>(DEFAULT_TABLE_SIZE);

		// the initial parameters
		try {
			setDefaultParams();
		} catch (PropertyVetoException ignored) {
			// because it never should happen in this context
			System.err.println(ignored.toString());
		}
	}

	/**
	 * Obtain an instance of this ...
	 * 
	 * @return atomInfo AtomProperty configuration
	 */
	public static AtomInfo getInstance() {
		if (atomInfo == null) {
			atomInfo = new AtomInfo();
		}

		return atomInfo;
	}

	/**
	 * method to reset the instance of AtomInfo .. so that default values are
	 * loaded. Note that the the values are not saved on to a persistent store.
	 * 
	 */
	public static void reset() {
		if (atomInfo == null)
			return;

		// re-initilize the parameters
		try {
			atomInfo.setDefaultParams();

			AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(atomInfo);

			changeEvent.setChangeType(AtomInfoChangeEvent.ALL_CANGED);
			// fire the event!
			atomInfo.fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
		} catch (PropertyVetoException ignored) {
			// because it never should happen in this context
			System.err.println(ignored.toString());
		} // end of try catch
	}

	/**
	 * Method to reset the instance of AtomInfo .. so that default values are loaded
	 * for a particular atom. Note that the values are not saved on to a persistent
	 * store.
	 * 
	 * @param symbol Atomic symbol
	 */
	public void resetValues(String symbol) {
		nameTable.put(symbol, originalNameTable.get(symbol));
		atomicNumberTable.put(symbol, originalAtomicNumberTable.get(symbol));
		atomicWeightTable.put(symbol, originalAtomicWeightTable.get(symbol));
		covalentRadiusTable.put(symbol, originalCovalentRadiusTable.get(symbol));
		vdwRadiusTable.put(symbol, originalVdwRadiusTable.get(symbol));
		defaultValencyTable.put(symbol, originalDefaultValencyTable.get(symbol));
		doubleBondOverlapTable.put(symbol, originalDoubleBondOverlapTable.get(symbol));
		weakBondAngleTable.put(symbol, originalWeakBondAngleTable.get(symbol));
		colorTable.put(symbol, originalColorTable.get(symbol));

		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ALL_CANGED);
		// fire the event!
		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Save user level configuration file
	 * 
	 * @throws IOException Signals that an I/O exception of some sort has occurred.
	 */
	public void saveUserConfigFile() throws IOException {
		// TODO: we are being a bit paranoid by replicating every thing, though
		// that is not really necessary
		StringResource strings = StringResource.getInstance();
		FileOutputStream fos = new FileOutputStream(strings.getUserAtomInfo());

		fos.write(strings.getXmlHeader().getBytes());

		fos.write("<atominfo> \n".getBytes());
		for (String key : nameTable.keySet()) {
			fos.write("\t<atom> \n".getBytes());
			fos.write(("\t\t<physical symbol=\"" + key + "\" name=\"" + nameTable.get(key) + "\" atomicNumber=\""
					+ atomicNumberTable.get(key) + "\" atomicWeight=\"" + atomicWeightTable.get(key)
					+ "\" covalentRadius=\"" + covalentRadiusTable.get(key) + "\" vdwRadius=\""
					+ vdwRadiusTable.get(key) + "\" defaultValency=\"" + defaultValencyTable.get(key)
					+ "\" weakBondAngle=\"" + weakBondAngleTable.get(key) + "\" doubleBondOverlap=\""
					+ doubleBondOverlapTable.get(key) + "\"/> \n").getBytes());
			fos.write("\t\t<display> \n".getBytes());
			Color theColor = colorTable.get(key);
			fos.write(("\t\t\t<color r=\"" + theColor.getRed() + "\" g=\"" + theColor.getGreen() + "\" b=\""
					+ theColor.getBlue() + "\"/>\n").getBytes());
			fos.write("\t\t</display>\n".getBytes());
			fos.write("\t</atom> \n".getBytes());
		} // end for
		fos.write("</atominfo> \n".getBytes());

		fos.close();
	}

	/**
	 * private method to set the default parameters
	 */
	private void setDefaultParams() throws PropertyVetoException {
		StringResource strings = StringResource.getInstance();
		try {
			// read the internal XML config file
			Document configDoc = Utility.parseXML(getClass().getResourceAsStream(strings.getDefaultAtomInfo()));

			// and save the info. properly
			saveOriginal(configDoc);

			// read in user XML config file
			// and override the settings
			saveUser();
		} catch (Exception e) {
			throw new PropertyVetoException("Exception in " + "AtomInfo.setParameter()" + e.toString(), null);
		} // end of try .. catch block
	}

	/**
	 * Recursive routine save DOM tree nodes
	 */
	private void saveOriginal(Node n) {
		int type = n.getNodeType(); // get node type

		switch (type) {
		case Node.ATTRIBUTE_NODE:
			String nodeName = n.getNodeName();

			if (nodeName.equals("name")) {
				originalNameTable.put(symbol, n.getNodeValue());
			} else if (nodeName.equals("atomicNumber")) {
				originalAtomicNumberTable.put(symbol, Integer.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("atomicWeight")) {
				originalAtomicWeightTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("covalentRadius")) {
				originalCovalentRadiusTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("vdwRadius")) {
				originalVdwRadiusTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("defaultValency")) {
				originalDefaultValencyTable.put(symbol, Integer.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("weakBondAngle")) {
				originalWeakBondAngleTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("doubleBondOverlap")) {
				originalDoubleBondOverlapTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} // end if

			break;
		case Node.ELEMENT_NODE:
			element = n.getNodeName();

			NamedNodeMap atts = n.getAttributes();
			if (element.equals("physical")) {
				symbol = atts.getNamedItem("symbol").getNodeValue();

				// save the others
				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveOriginal(att);
				} // end for
			} else if (element.equals("color")) {
				originalColorTable.put(symbol,
						new Color(Integer.parseInt(atts.getNamedItem("r").getNodeValue()),
								Integer.parseInt(atts.getNamedItem("g").getNodeValue()),
								Integer.parseInt(atts.getNamedItem("b").getNodeValue())));
			} else {
				if (atts == null)
					return;

				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveOriginal(att);
				} // end for
			} // end if

			break;
		default:
			break;
		} // end switch..case

		// save children if any
		for (Node child = n.getFirstChild(); child != null; child = child.getNextSibling()) {
			saveOriginal(child);
		} // end for
	} // end of method saveOriginal()

	/**
	 * Recursive routine save DOM tree nodes
	 */
	private void saveUserNode(Node n) {
		int type = n.getNodeType(); // get node type

		switch (type) {
		case Node.ATTRIBUTE_NODE:
			String nodeName = n.getNodeName();

			if (nodeName.equals("name")) {
				nameTable.put(symbol, n.getNodeValue());
			} else if (nodeName.equals("atomicNumber")) {
				atomicNumberTable.put(symbol, Integer.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("atomicWeight")) {
				atomicWeightTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("covalentRadius")) {
				covalentRadiusTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("vdwRadius")) {
				vdwRadiusTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("defaultValency")) {
				defaultValencyTable.put(symbol, Integer.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("weakBondAngle")) {
				weakBondAngleTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} else if (nodeName.equals("doubleBondOverlap")) {
				doubleBondOverlapTable.put(symbol, Double.valueOf(n.getNodeValue()));
			} // end if

			break;
		case Node.ELEMENT_NODE:
			element = n.getNodeName();

			NamedNodeMap atts = n.getAttributes();
			if (element.equals("physical")) {
				symbol = atts.getNamedItem("symbol").getNodeValue();

				// save the others
				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveUserNode(att);
				} // end for
			} else if (element.equals("color")) {
				colorTable.put(symbol,
						new Color(Integer.parseInt(atts.getNamedItem("r").getNodeValue()),
								Integer.parseInt(atts.getNamedItem("g").getNodeValue()),
								Integer.parseInt(atts.getNamedItem("b").getNodeValue())));
			} else {
				if (atts == null)
					return;

				for (int i = 0; i < atts.getLength(); i++) {
					Node att = atts.item(i);
					saveUserNode(att);
				} // end for
			} // end if

			break;
		default:
			break;
		} // end switch..case

		// save children if any
		for (Node child = n.getFirstChild(); child != null; child = child.getNextSibling()) {
			saveUserNode(child);
		} // end for
	} // end of method saveUserNode()

	/** make a copy of the original */
	private void copyOriginals() {
		// just make a copy of original
		copyTableS(originalNameTable, nameTable);
		copyTableI(originalAtomicNumberTable, atomicNumberTable);
		copyTableD(originalAtomicWeightTable, atomicWeightTable);
		copyTableD(originalCovalentRadiusTable, covalentRadiusTable);
		copyTableD(originalVdwRadiusTable, vdwRadiusTable);
		copyTableI(originalDefaultValencyTable, defaultValencyTable);
		copyTableD(originalDoubleBondOverlapTable, doubleBondOverlapTable);
		copyTableD(originalWeakBondAngleTable, weakBondAngleTable);
		copyTableC(originalColorTable, colorTable);
	}

	/**
	 * Routine save DOM tree nodes
	 */
	private void saveUser() throws Exception {
		StringResource strings = StringResource.getInstance();

		try {
			if ((new File(strings.getUserAtomInfo())).exists()) {
				// if a user configuration file exists, then we copy
				// its contents
				Document configDoc = Utility.parseXML(strings.getUserAtomInfo());

				saveUserNode(configDoc);
			} else {
				copyOriginals();
			} // end if
		} catch (Exception e) {
			System.err.println("Restricted mode copying.");
			copyOriginals();
		} // end if
	} // end of method saveUser()

	/**
	 * Copy one table to another (of strings)
	 */
	private void copyTableS(HashMap<String, String> src, HashMap<String, String> dest) {
		for (String key : src.keySet()) {
			dest.put(key, src.get(key));
		} // end for
	}

	/**
	 * Copy one table to another (of ints)
	 */
	private void copyTableI(HashMap<String, Integer> src, HashMap<String, Integer> dest) {
		for (String key : src.keySet()) {
			dest.put(key, src.get(key));
		} // end for
	}

	/**
	 * Copy one table to another (of doubles)
	 */
	private void copyTableD(HashMap<String, Double> src, HashMap<String, Double> dest) {
		for (String key : src.keySet()) {
			dest.put(key, src.get(key));
		} // end for
	}

	/**
	 * Copy one table to another (of Color)
	 */
	private void copyTableC(HashMap<String, Color> src, HashMap<String, Color> dest) {
		for (String key : src.keySet()) {
			dest.put(key, src.get(key));
		} // end for
	}

	/**
	 * Method returns the value of parameter pertaining to the key.
	 * 
	 * @return Parameter - the parameter value.
	 * @throws NullPointerException if the key is not found
	 */
	@Override
	public Parameter getParameter(String key) {
		return new AtomProperty((String) nameTable.get(key), (atomicNumberTable.get(key)).intValue(),
				(atomicWeightTable.get(key)).doubleValue(), (defaultValencyTable.get(key)).intValue(),
				(covalentRadiusTable.get(key)).doubleValue(), (vdwRadiusTable.get(key)).doubleValue(),
				(weakBondAngleTable.get(key)).doubleValue(), (doubleBondOverlapTable.get(key)).doubleValue(),
				(Color) colorTable.get(key));
	}

	/**
	 * method to set the new parameter value.
	 * 
	 * @param key       - the key whose value needs to be changed or added
	 * @param parameter - the new parameter value
	 * @throws PropertyVetoException - incase changing property value is not
	 *                               supported
	 */
	@Override
	public void setParameter(String key, Parameter parameter) throws PropertyVetoException {
		AtomProperty ap = (AtomProperty) parameter;

		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ALL_CANGED);

		changeEvent.setOldValue(new AtomProperty(Utility.getString(nameTable, key),
				Utility.getInteger(atomicNumberTable, key), Utility.getDouble(atomicWeightTable, key),
				Utility.getInteger(defaultValencyTable, key), Utility.getDouble(covalentRadiusTable, key),
				Utility.getDouble(vdwRadiusTable, key), Utility.getDouble(weakBondAngleTable, key),
				Utility.getDouble(doubleBondOverlapTable, key), Utility.getColor(colorTable, key))); // old val
		changeEvent.setNewValue(ap); // new val
		changeEvent.setAtomSymbol(key); // the atomic symbol

		// change to new value
		nameTable.put(key, ap.getName());
		atomicNumberTable.put(key, Integer.valueOf(ap.getAtomicNumber()));
		atomicWeightTable.put(key, Double.valueOf(ap.getAtomicWeight()));
		defaultValencyTable.put(key, Integer.valueOf(ap.getDefaultValency()));
		covalentRadiusTable.put(key, Double.valueOf(ap.getCovalentRadius()));
		vdwRadiusTable.put(key, Double.valueOf(ap.getVdwRadius()));
		weakBondAngleTable.put(key, Double.valueOf(ap.getWeakBondAngle()));
		doubleBondOverlapTable.put(key, Double.valueOf(ap.getDoubleBondOverlap()));
		colorTable.put(key, ap.getColor());

		// fire the event!
		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * is the configuration stored as file?
	 * 
	 * @return boolean - true / false
	 */
	@Override
	public boolean isStoredAsFile() {
		return true;
	}

	/**
	 * Method returns the absolute path of the configuration file if
	 * isStoredAsFile() call results in a true value, else a null is returned.
	 * 
	 * @return String - the absolute path of the configuration file
	 */
	@Override
	public String getConfigFile() {
		return StringResource.getInstance().getDefaultAtomInfo();
	}

	/**
	 * changes current value of file storage
	 * 
	 * @param storedAsFile new value of this storage
	 */
	@Override
	public void setStoredAsFile(boolean storedAsFile) {
		// igonored! ... we don't want to listen here
	}

	/**
	 * sets the value of a parameter
	 * 
	 * @param configFile - the new name of the config file
	 * @throws PropertyVetoException some problem can't change the stuff!
	 */
	@Override
	public void setConfigFile(String configFile) throws PropertyVetoException {
		throw new PropertyVetoException("Cannot change this property!", null);
	}

	/**
	 * Getter for property nameTable.
	 * 
	 * @return Value of property nameTable.
	 * 
	 */
	public HashMap<String, String> getNameTable() {
		return this.nameTable;
	}

	/**
	 * Setter for property nameTable.
	 * 
	 * @param nameTable New value of property nameTable.
	 * 
	 */
	public void setNameTable(HashMap<String, String> nameTable) {
		this.nameTable = nameTable;
	}

	/**
	 * Getter for property atomicNumberTable.
	 * 
	 * @return Value of property atomicNumberTable.
	 * 
	 */
	public HashMap<String, Integer> getAtomicNumberTable() {
		return this.atomicNumberTable;
	}

	/**
	 * Setter for property atomicNumberTable.
	 * 
	 * @param atomicNumberTable New value of property atomicNumberTable.
	 * 
	 */
	public void setAtomicNumberTable(HashMap<String, Integer> atomicNumberTable) {
		this.atomicNumberTable = atomicNumberTable;
	}

	/**
	 * Getter for property atomicWeightTable.
	 * 
	 * @return Value of property atomicWeightTable.
	 * 
	 */
	public HashMap<String, Double> getAtomicWeightTable() {
		return this.atomicWeightTable;
	}

	/**
	 * Setter for property atomicWeightTable.
	 * 
	 * @param atomicWeightTable New value of property atomicWeightTable.
	 * 
	 */
	public void setAtomicWeightTable(HashMap<String, Double> atomicWeightTable) {
		this.atomicWeightTable = atomicWeightTable;
	}

	/**
	 * Getter for property covalentRadiusTable.
	 * 
	 * @return Value of property covalentRadiusTable.
	 * 
	 */
	public HashMap<String, Double> getCovalentRadiusTable() {
		return this.covalentRadiusTable;
	}

	/**
	 * Setter for property covalentRadiusTable.
	 * 
	 * @param covalentRadiusTable New value of property covalentRadiusTable.
	 * 
	 */
	public void setCovalentRadiusTable(HashMap<String, Double> covalentRadiusTable) {
		this.covalentRadiusTable = covalentRadiusTable;
	}

	/**
	 * Getter for property vdwRadiusTable.
	 * 
	 * @return Value of property vdwRadiusTable.
	 * 
	 */
	public HashMap<String, Double> getVdwRadiusTable() {
		return this.vdwRadiusTable;
	}

	/**
	 * Setter for property vdwRadiusTable.
	 * 
	 * @param vdwRadiusTable New value of property vdwRadiusTable.
	 * 
	 */
	public void setVdwRadiusTable(HashMap<String, Double> vdwRadiusTable) {
		this.vdwRadiusTable = vdwRadiusTable;
	}

	/**
	 * Getter for property defaultValencyTable.
	 * 
	 * @return Value of property defaultValencyTable.
	 * 
	 */
	public HashMap<String, Integer> getDefaultValencyTable() {
		return this.defaultValencyTable;
	}

	/**
	 * Setter for property defaultValencyTable.
	 * 
	 * @param defaultValencyTable New value of property defaultValencyTable.
	 * 
	 */
	public void setDefaultValencyTable(HashMap<String, Integer> defaultValencyTable) {
		this.defaultValencyTable = defaultValencyTable;
	}

	/**
	 * Getter for property weakBondAngleTable.
	 * 
	 * @return Value of property weakBondAngleTable.
	 * 
	 */
	public HashMap<String, Double> getWeakBondAngleTable() {
		return this.weakBondAngleTable;
	}

	/**
	 * Setter for property weakBondAngleTable.
	 * 
	 * @param weakBondAngleTable New value of property weakBondAngleTable.
	 * 
	 */
	public void setWeakBondAngleTable(HashMap<String, Double> weakBondAngleTable) {
		this.weakBondAngleTable = weakBondAngleTable;
	}

	/**
	 * Getter for property doubleBondOverlap.
	 * 
	 * @return Value of property doubleBondOverlap.
	 * 
	 */
	public HashMap<String, Double> getDoubleBondOverlapTable() {
		return this.doubleBondOverlapTable;
	}

	/**
	 * Setter for property doubleBondOverlap.
	 * 
	 * @param doubleBondOverlapTable New value of property doubleBondOverlapTable.
	 * 
	 */
	public void setDoubleBondOverlapTable(HashMap<String, Double> doubleBondOverlapTable) {
		this.doubleBondOverlapTable = doubleBondOverlapTable;
	}

	/**
	 * Getter for property colorTable.
	 * 
	 * @return Value of property colorTable.
	 * 
	 */
	public HashMap<String, Color> getColorTable() {
		return this.colorTable;
	}

	/**
	 * Setter for property colorTable.
	 * 
	 * @param colorTable New value of property colorTable.
	 * 
	 */
	public void setColorTable(HashMap<String, Color> colorTable) {
		this.colorTable = colorTable;
	}

	/**
	 * Registers AtomInfoChangeListener to receive events.
	 * 
	 * @param listener The listener to register.
	 * 
	 */
	public synchronized void addAtomInfoChangeListener(AtomInfoChangeListener listener) {
		if (listenerList == null) {
			listenerList = new EventListenerList<>();
		}
		listenerList.add(AtomInfoChangeListener.class, listener);
	}

	/**
	 * Removes AtomInfoChangeListener from the list of listeners.
	 * 
	 * @param listener The listener to remove.
	 * 
	 */
	public synchronized void removeAtomInfoChangeListener(AtomInfoChangeListener listener) {
		listenerList.remove(AtomInfoChangeListener.class, listener);
	}

	/**
	 * Notifies all registered listeners about the event.
	 * 
	 * @param event The event to be fired
	 * 
	 */
	private void fireAtomInfoChangeListenerAtomInfoChanged(AtomInfoChangeEvent event) {
		if (listenerList == null)
			return;

		// the great hack ... bypass notification if there was no real
		// change!
		if (event.getOldValue() != null) {
			if (event.getOldValue().equals(event.getNewValue()))
				return;
		} // end if

		for (Object listener : listenerList.getListenerList()) {
			((AtomInfoChangeListener) listener).atomInfoChanged(event);
		} // end for

	}

	/**
	 * Getter for property atomicNumber.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property atomicNumber for the specified symbol
	 */
	public int getAtomicNumber(String symbol) {
		try {
			return atomicNumberTable.get(symbol);
		} catch (Exception e) {
			return atomicNumberTable.get("X");
		} // end of try catch block
	}

	/**
	 * Getter for property symbol.
	 * 
	 * @param atomicNumber for the required symbol
	 * @return symbol - the atom symbol, IUPAC name! If no such atomic number is
	 *         found then "X" is returned.
	 */
	public String getSymbol(int atomicNumber) {
		Enumeration<Integer> eles = (Enumeration<Integer>) atomicNumberTable.values();
		Enumeration<String> keys = (Enumeration<String>) atomicNumberTable.keySet();

		while (eles.hasMoreElements()) {
			if (eles.nextElement() == atomicNumber)
				return keys.nextElement();
			keys.nextElement();
		} // end while

		return "X";
	}

	/**
	 * Setter for property atomicNumber.
	 * 
	 * @param symbol       - the atom symbol, IUPAC name!
	 * @param atomicNumber New value of property atomicNumber.
	 */
	public void setAtomicNumber(String symbol, int atomicNumber) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ATOMIC_NUMBER);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(atomicNumberTable.get(symbol));

		atomicNumberTable.put(symbol, Integer.valueOf(atomicNumber));

		changeEvent.setNewValue(Integer.valueOf(atomicNumber));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property atomicWeight.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property atomicWeight for the specified symbol
	 */
	public double getAtomicWeight(String symbol) {
		try {
			return atomicWeightTable.get(symbol);
		} catch (Exception e) {
			return atomicWeightTable.get("X");
		} // end of try .. catch block
	}

	/**
	 * Setter for property atomicWeight.
	 * 
	 * @param symbol       - the atom symbol, IUPAC name!
	 * @param atomicWeight New value of property atomicWeight.
	 */
	public void setAtomicWeight(String symbol, double atomicWeight) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ATOMIC_WEIGHT);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(atomicWeightTable.get(symbol));

		atomicWeightTable.put(symbol, Double.valueOf(atomicWeight));

		changeEvent.setNewValue(Double.valueOf(atomicWeight));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property defaultValency.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property defaultValency for the specified symbol
	 */
	public int getDefaultValency(String symbol) {
		try {
			return defaultValencyTable.get(symbol);
		} catch (Exception e) {
			return defaultValencyTable.get("X");
		} // end of try .. catch block
	}

	/**
	 * Setter for property defaultValency.
	 * 
	 * @param symbol         - the atom symbol, IUPAC name!
	 * @param defaultValency New value of property defaultValency.
	 */
	public void setDefaultValency(String symbol, int defaultValency) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.DEFAULT_VALENCY);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(defaultValencyTable.get(symbol));

		defaultValencyTable.put(symbol, Integer.valueOf(defaultValency));

		changeEvent.setNewValue(Integer.valueOf(defaultValency));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property covalentRadius.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property covalentRadius for the specified symbol
	 */
	public double getCovalentRadius(String symbol) {
		try {
			return covalentRadiusTable.get(symbol);
		} catch (Exception e) {
			return covalentRadiusTable.get("X");
		} // end of try ... catch block
	}

	/**
	 * Setter for property covalentRadius.
	 * 
	 * @param symbol         - the atom symbol, IUPAC name!
	 * @param covalentRadius New value of property covalentRadius.
	 */
	public void setCovalentRadius(String symbol, double covalentRadius) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.COVALENT_RADIUS);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(covalentRadiusTable.get(symbol));

		covalentRadiusTable.put(symbol, Double.valueOf(covalentRadius));

		changeEvent.setNewValue(Double.valueOf(covalentRadius));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property vdwRadius.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property vdwRadius for the specified symbol
	 */
	public double getVdwRadius(String symbol) {
		try {
			return vdwRadiusTable.get(symbol);
		} catch (Exception e) {
			return vdwRadiusTable.get("X");
		} // end try .. catch block
	}

	/**
	 * Setter for property vdwRadius.
	 * 
	 * @param symbol    - the atom symbol, IUPAC name!
	 * @param vdwRadius New value of property vdwRadius.
	 */
	public void setVdwRadius(String symbol, double vdwRadius) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.VDW_RADIUS);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(vdwRadiusTable.get(symbol));

		vdwRadiusTable.put(symbol, Double.valueOf(vdwRadius));

		changeEvent.setNewValue(Double.valueOf(vdwRadius));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property weakBondAngle.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property weakBondAngle for the specified symbol
	 */
	public double getWeakBondAngle(String symbol) {
		try {
			return weakBondAngleTable.get(symbol);
		} catch (Exception e) {
			return weakBondAngleTable.get("X");
		} // end of try .. catch block
	}

	/**
	 * Setter for property weakBondAngle.
	 * 
	 * @param symbol        - the atom symbol, IUPAC name!
	 * @param weakBondAngle New value of property weakBondAngle.
	 */
	public void setWeakBondAngle(String symbol, double weakBondAngle) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.WEAK_BOND_ANGLE);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(weakBondAngleTable.get(symbol));

		weakBondAngleTable.put(symbol, Double.valueOf(weakBondAngle));

		changeEvent.setNewValue(Double.valueOf(weakBondAngle));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property name.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property name for the specified symbol
	 */
	public String getName(String symbol) {
		try {
			return nameTable.get(symbol);
		} catch (Exception e) {
			return nameTable.get("X");
		} // end of try ... catch block
	}

	/**
	 * Setter for property weakBondAngle.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @param name   New value of property name.
	 */
	public void setName(String symbol, String name) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.ATOM_NAME);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(nameTable.get(symbol));

		nameTable.put(symbol, name);

		changeEvent.setNewValue(name);

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property color.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property color for the specified symbol
	 */
	public Color getColor(String symbol) {
		try {
			Color clr = colorTable.get(symbol);

			if (clr == null)
				return colorTable.get("X");

			return clr;
		} catch (Exception e) {
			return colorTable.get("X");
		} // end of try .. catch block
	}

	/**
	 * Setter for property color.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @param color  New value of property color.
	 */
	public void setColor(String symbol, Color color) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.COLOR_VALUE);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(colorTable.get(symbol));

		colorTable.put(symbol, color);

		changeEvent.setNewValue(color);

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}

	/**
	 * Getter for property doubleBondOverlap.
	 * 
	 * @param symbol - the atom symbol, IUPAC name!
	 * @return Value of property doubleBondOverlap for the specified symbol
	 */
	public double getDoubleBondOverlap(String symbol) {
		try {
			return doubleBondOverlapTable.get(symbol);
		} catch (Exception e) {
			return doubleBondOverlapTable.get("X");
		} // end of try ... catch block
	}

	/**
	 * Setter for property doubleBondOverlap.
	 * 
	 * @param symbol            - the atom symbol, IUPAC name!
	 * @param doubleBondOverlap New value of property doubleBondOverlap.
	 */
	public void setDoubleBondOverlap(String symbol, double doubleBondOverlap) {
		AtomInfoChangeEvent changeEvent = new AtomInfoChangeEvent(this);

		changeEvent.setChangeType(AtomInfoChangeEvent.DOUBLE_BOND_OVERLAP);
		changeEvent.setAtomSymbol(symbol);
		changeEvent.setOldValue(doubleBondOverlapTable.get(symbol));

		doubleBondOverlapTable.put(symbol, Double.valueOf(doubleBondOverlap));

		changeEvent.setNewValue(Double.valueOf(doubleBondOverlap));

		fireAtomInfoChangeListenerAtomInfoChanged(changeEvent);
	}
}
