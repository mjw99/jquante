package name.mjw.jquante.math.qm.basis;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.Iterator;

import org.junit.jupiter.api.Test;

class PowerListTest {

    @Test
    void getInstanceReturnsSingleton() {
        PowerList p1 = PowerList.getInstance();
        PowerList p2 = PowerList.getInstance();
        assertNotNull(p1);
        assertSame(p1, p2);
    }

    @Test
    void generatePowerListSHasOneElement() {
        ArrayList<Power> list = PowerList.getInstance().generatePowerList(0);
        assertEquals(1, list.size());
        assertEquals(new Power(0, 0, 0), list.get(0));
    }

    @Test
    void generatePowerListPHasThreeElementsWithCorrectAngularMomentum() {
        ArrayList<Power> list = PowerList.getInstance().generatePowerList(1);
        assertEquals(3, list.size());
        for (Power p : list) {
            assertEquals(1, p.getTotalAngularMomentum());
        }
    }

    @Test
    void generatePowerListDHasSixElements() {
        ArrayList<Power> list = PowerList.getInstance().generatePowerList(2);
        assertEquals(6, list.size());
        for (Power p : list) {
            assertEquals(2, p.getTotalAngularMomentum());
        }
    }

    @Test
    void generatePowerListFHasTenElements() {
        ArrayList<Power> list = PowerList.getInstance().generatePowerList(3);
        assertEquals(10, list.size());
        for (Power p : list) {
            assertEquals(3, p.getTotalAngularMomentum());
        }
    }

    @Test
    void generatePowerListGHasFifteenElements() {
        ArrayList<Power> list = PowerList.getInstance().generatePowerList(4);
        assertEquals(15, list.size());
    }

    @Test
    void generatePowerListHHasTwentyOneElements() {
        ArrayList<Power> list = PowerList.getInstance().generatePowerList(5);
        assertEquals(21, list.size());
    }

    @Test
    void getPowerListSIteratorHasOneElement() {
        Iterator<Power> iter = PowerList.getInstance().getPowerList("S");
        assertTrue(iter.hasNext());
        Power p = iter.next();
        assertEquals(0, p.getTotalAngularMomentum());
        assertFalse(iter.hasNext());
    }

    @Test
    void getPowerListPIteratorHasThreeElements() {
        Iterator<Power> iter = PowerList.getInstance().getPowerList("P");
        int count = 0;
        while (iter.hasNext()) {
            iter.next();
            count++;
        }
        assertEquals(3, count);
    }

    @Test
    void getPowerListDIteratorHasSixElements() {
        Iterator<Power> iter = PowerList.getInstance().getPowerList("D");
        int count = 0;
        while (iter.hasNext()) {
            iter.next();
            count++;
        }
        assertEquals(6, count);
    }
}
