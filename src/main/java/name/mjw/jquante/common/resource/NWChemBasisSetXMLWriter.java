package name.mjw.jquante.common.resource;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

/**
 * 
 * @author mw529
 * @see <a href="http://www.vogella.com/tutorials/JAXB/article.html">Example
 *      tutorial</a>
 */
public class NWChemBasisSetXMLWriter {

	public static void main(String[] args) throws JAXBException {

		String fileName = "/usr/local/shared/ubuntu-16.04/x86_64/nwchem/6.6/data/libraries/aug-cc-pvtz";

		NWChemBasisSetFile nWChemBasisSetFile = new NWChemBasisSetFile();

		nWChemBasisSetFile.read(fileName);

		// create JAXB context and instantiate marshaller
		JAXBContext context = null;
		context = JAXBContext.newInstance(NWChemBasisSetFile.class);

		Marshaller m = context.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

		// Write to System.out
		m.marshal(nWChemBasisSetFile, System.out);

	}

}
