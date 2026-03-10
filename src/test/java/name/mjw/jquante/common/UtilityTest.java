package name.mjw.jquante.common;

import static org.junit.jupiter.api.Assertions.*;

import java.io.ByteArrayInputStream;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;

import org.junit.jupiter.api.Test;
import org.w3c.dom.Document;

class UtilityTest {

    @Test
    void capitaliseFirstLetterOnly() {
        assertEquals("Carbon", Utility.capitalise("carbon"));
    }

    @Test
    void capitaliseAlreadyCapitalised() {
        assertEquals("Hydrogen", Utility.capitalise("Hydrogen"));
    }

    @Test
    void capitaliseSingleCharacter() {
        assertEquals("H", Utility.capitalise("h"));
    }

    @Test
    void capitalisePreservesRestOfString() {
        assertEquals("He", Utility.capitalise("he"));
    }

    @Test
    void getIntegerReturnsValueWhenPresent() {
        HashMap<String, Integer> map = new HashMap<>();
        map.put("key", 42);
        assertEquals(42, Utility.getInteger(map, "key"));
    }

    @Test
    void getIntegerReturnsZeroWhenMissing() {
        HashMap<String, Integer> map = new HashMap<>();
        assertEquals(0, Utility.getInteger(map, "missing"));
    }

    @Test
    void getDoubleReturnsValueWhenPresent() {
        HashMap<String, Double> map = new HashMap<>();
        map.put("pi", 3.14);
        assertEquals(3.14, Utility.getDouble(map, "pi"), 1e-10);
    }

    @Test
    void getDoubleReturnsZeroWhenMissing() {
        HashMap<String, Double> map = new HashMap<>();
        assertEquals(0.0, Utility.getDouble(map, "missing"), 1e-10);
    }

    @Test
    void getStringReturnsValueWhenPresent() {
        HashMap<String, String> map = new HashMap<>();
        map.put("greeting", "hello");
        assertEquals("hello", Utility.getString(map, "greeting"));
    }

    @Test
    void getStringReturnsBlankWhenMissing() {
        HashMap<String, String> map = new HashMap<>();
        assertEquals(" ", Utility.getString(map, "missing"));
    }

    @Test
    void parseXMLReturnsDocument() throws Exception {
        String xml = "<root><child/></root>";
        ByteArrayInputStream stream = new ByteArrayInputStream(xml.getBytes(StandardCharsets.UTF_8));
        Document doc = Utility.parseXML(stream);
        assertNotNull(doc);
        assertEquals("root", doc.getDocumentElement().getTagName());
    }

    @Test
    void auToAngstromFactorIsPositive() {
        assertTrue(Utility.AU_TO_ANGSTROM_FACTOR > 0.0);
    }

    @Test
    void kcalToAuFactorIsPositive() {
        assertTrue(Utility.KCAL_TO_AU_FACTOR > 0.0);
    }
}
