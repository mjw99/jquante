package name.mjw.jquante.common.resource;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

/**
 * Command-line utility that reads a NWChem basis-set library file and writes it
 * out as an XML file using JAXB marshalling.
 *
 * @author mw529
 * @see <a href="http://www.vogella.com/tutorials/JAXB/article.html">Example JAXB tutorial</a>
 */
public class NWChemBasisSetXMLWriter {

	/**
	 * Entry point. Reads the NWChem basis file at the hard-coded path, marshals
	 * it to XML, and writes it to {@code basis_<name>.xml} in the current
	 * working directory.
	 *
	 * @param args command-line arguments (not used).
	 * @throws JAXBException if JAXB context creation or marshalling fails.
	 */
	public static void main(String[] args) throws JAXBException {


		String fileName = "/usr/local/shared/ubuntu-16.04/x86_64/nwchem/6.6/data/libraries/sto-3g";
		
		NWChemBasisSetFile nWChemBasisSetFile = new NWChemBasisSetFile();

		nWChemBasisSetFile.read(fileName);

		// create JAXB context and instantiate marshaller
		JAXBContext context = null;
		context = JAXBContext.newInstance(NWChemBasisSetFile.class);

		Marshaller m = context.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

		try {
			OutputStream os = new FileOutputStream(
					"basis_" + nWChemBasisSetFile.getBasisSetName() + ".xml");
			m.marshal(nWChemBasisSetFile, os);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

}
