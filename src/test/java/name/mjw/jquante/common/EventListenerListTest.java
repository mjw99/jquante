package name.mjw.jquante.common;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class EventListenerListTest {

    interface SimpleListener {
        void onEvent();
    }

    @SuppressWarnings("unchecked")
    private EventListenerList<SimpleListener> list;

    @BeforeEach
    @SuppressWarnings("unchecked")
    void setUp() {
        list = new EventListenerList<>();
    }

    @Test
    void emptyListHasZeroListeners() {
        assertEquals(0, list.getListenerCount());
    }

    @Test
    void addOneListenerIncreasesCount() {
        list.add(SimpleListener.class, () -> {});
        assertEquals(1, list.getListenerCount());
    }

    @Test
    void addTwoListenersCountIsTwo() {
        list.add(SimpleListener.class, () -> {});
        list.add(SimpleListener.class, () -> {});
        assertEquals(2, list.getListenerCount());
    }

    @Test
    void removeListenerDecreasesCount() {
        SimpleListener l = () -> {};
        list.add(SimpleListener.class, l);
        list.remove(SimpleListener.class, l);
        assertEquals(0, list.getListenerCount());
    }

    @Test
    void removeFromEmptyListDoesNotThrow() {
        assertDoesNotThrow(() -> list.remove(SimpleListener.class, () -> {}));
    }

    @Test
    void getListenerListReturnsAllListeners() {
        SimpleListener l1 = () -> {};
        SimpleListener l2 = () -> {};
        list.add(SimpleListener.class, l1);
        list.add(SimpleListener.class, l2);
        assertEquals(2, list.getListenerList().length);
    }

    @Test
    void getListenerListByClassReturnsMatchingListeners() {
        SimpleListener l = () -> {};
        list.add(SimpleListener.class, l);
        Object[] result = list.getListenerList(SimpleListener.class);
        assertEquals(1, result.length);
        assertSame(l, result[0]);
    }

    @Test
    void getListenerCountByClassReturnsCorrectCount() {
        list.add(SimpleListener.class, () -> {});
        list.add(SimpleListener.class, () -> {});
        assertEquals(2, list.getListenerCount(SimpleListener.class));
    }

    @Test
    void getListenerCountByUnregisteredClassReturnsZero() {
        assertEquals(0, list.getListenerCount(SimpleListener.class));
    }

    @Test
    void getListenerListByUnregisteredClassReturnsEmpty() {
        assertEquals(0, list.getListenerList(SimpleListener.class).length);
    }
}
