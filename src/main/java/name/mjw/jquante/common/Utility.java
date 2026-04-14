package name.mjw.jquante.common;

import org.w3c.dom.Document;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

/**
 * Utility class providing basic functionality such as ID generation etc. The
 * class it self is declared as final, and all the methods exposed are static.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class Utility {

	/**
	 * 1 a.u. = 0.52917724924 angstroms When converting from a.u. to angstroms,
	 * multiply by this factor. When converting from angstroms to a.u., divide by
	 * this factor.
	 */
	public static final double AU_TO_ANGSTROM_FACTOR = 0.52917724924;

	/**
	 * 1 kcal/mol = 0.00159362 AU.
	 */
	public static final double KCAL_TO_AU_FACTOR = 0.00159362;

	/**
	 * private constructor ... avoid object creation
	 */
	private Utility() {
	}

	/**
	 * Method to capitalise a string.
	 *
	 * @param theString the string to be capitalised
	 * @return the capitalised string
	 * @throws NullPointerException      if {@code theString} is null.
	 * @throws StringIndexOutOfBoundsException if {@code theString} is empty.
	 */
	public static String capitalise(String theString) {
		String firstChar = theString.substring(0, 1);
		String restString = theString.substring(1, theString.length());

		if (restString == null)
			return firstChar.toUpperCase();
		else
			return (firstChar.toUpperCase() + restString);
	}

	/**
	 * A simple utility method to parse XML using a default DOM parser and return
	 * the results as an abstract instance of Document class.
	 *
	 * @param xmlFile fully qualified path of the XML file to be parsed.
	 * @return the parsed XML document.
	 * @throws ParserConfigurationException if the parser cannot be configured.
	 * @throws SAXException                 if a parse error occurs.
	 * @throws IOException                  if the file cannot be read or is not found.
	 */
	public static Document parseXML(String xmlFile) throws ParserConfigurationException, SAXException, IOException {
		return parseXML(new FileInputStream(xmlFile));
	}

	/**
	 * A simple utility method to parse XML using a default DOM parser and return
	 * the results as an abstract instance of Document class.
	 *
	 * @param xmlStream an input stream supplying the XML data.
	 * @return the parsed XML document.
	 * @throws ParserConfigurationException if the parser cannot be configured.
	 * @throws SAXException                 if a parse error occurs.
	 * @throws IOException                  if the stream cannot be read.
	 */
	public static Document parseXML(InputStream xmlStream)
			throws ParserConfigurationException, SAXException, IOException {
		// we use the JAXP DOM parser for the job!
		// get an instance of the parser
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();

		// set the error handler .. to detect errors(?) in xml
		// in normal situations, error should never occur because
		// the xml files parsed here are automatically generated ones
		db.setErrorHandler(new ErrorHandler() {
			@Override
			public void warning(SAXParseException e) throws SAXException {
				System.err.println("Warning : " + e);
				e.printStackTrace();
			}

			@Override
			public void error(SAXParseException e) throws SAXException {
				System.err.println("Error : " + e);
				e.printStackTrace();
			}

			@Override
			public void fatalError(SAXParseException e) throws SAXException {
				System.err.println("Error : " + e);
				e.printStackTrace();
			}
		}); // setErrorHandler

		// and now, parse the input file:
		return db.parse(xmlStream);
	}

	/**
	 * Returns an integer value from a map. If the key is not present, zero is
	 * returned.
	 *
	 * @param table the map containing Integer values
	 * @param key   the key to look up
	 * @return the int value for the key, or zero if the key is not found
	 */
	public static int getInteger(HashMap<?, ?> table, Object key) {
		int value = 0;

		if (table.containsKey(key)) {
			value = ((Integer) table.get(key)).intValue();
		}

		return value;
	}

	/**
	 * Returns a double value from a map. If the key is not present, zero is
	 * returned.
	 *
	 * @param table the map containing Double values
	 * @param key   the key to look up
	 * @return the double value for the key, or zero if the key is not found
	 */
	public static double getDouble(HashMap<?, ?> table, Object key) {
		double value = 0.0;

		if (table.containsKey(key)) {
			value = ((Double) table.get(key)).doubleValue();
		}
		return value;
	}

	/**
	 * Returns a String value from a map. If the key is not present, a blank string
	 * is returned.
	 *
	 * @param table the map containing String values
	 * @param key   the key to look up
	 * @return the String value for the key, or a blank string if the key is not found
	 */
	public static String getString(HashMap<?, ?> table, Object key) {
		String value = " ";

		if (table.containsKey(key)) {
			value = (String) table.get(key);
		}
		return value;
	}
}
