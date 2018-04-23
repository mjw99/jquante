package name.mjw.jquante.parallel;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A simple framework for execution of tasks in parallel on a shared memory
 * system.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class SimpleParallelTaskExecuter {

	protected int noOfProcessors;

	private ArrayList<SimpleParallelTask> pTaskList;

	/** Creates a new instance of SimpleParallelTaskExecuter */
	public SimpleParallelTaskExecuter() {
		noOfProcessors = Runtime.getRuntime().availableProcessors()*4;
	}

	/**
	 * Execute the parallel task depending on the number of CPUs available
	 * 
	 * @param pTask
	 *            Simple parallel task.
	 */
	public void execute(SimpleParallelTask pTask) {
		ArrayList<Thread> parallelTasks;
		// init lists
		parallelTasks = new ArrayList<>(noOfProcessors);
		pTaskList = new ArrayList<>(noOfProcessors);

		// init items to distribute
		int totalItems = pTask.getTotalItems();
		int[] itemsPerProcessor = new int[noOfProcessors];

		while (totalItems > 0) {
			for (int j = 0; j < noOfProcessors; j++) {
				itemsPerProcessor[j] += (totalItems <= 0) ? 0 : 1;
				totalItems--;
				if (totalItems < 0)
					break;
			}
		}

		// create objects representing the distribution
		int i;
		int curIndx = 0;
		int start;
		int end;
		for (i = 0; i < noOfProcessors; i++) {
			start = curIndx;
			end = start + itemsPerProcessor[i];
			pTaskList.add(pTask.init(start, end));
			curIndx = end;
		}

		// then package them as threads, and start executing them
		for (i = 0; i < noOfProcessors; i++) {
			SimpleParallelTask spTask = pTaskList.get(i);
			spTask.setTotalItems(pTask.getTotalItems());
			spTask.setTaskName(pTask.getTaskName());

			Thread pTaskThread = new Thread(spTask);
			pTaskThread.setName(spTask.getTaskName());
			pTaskThread.start();

			parallelTasks.add(pTaskThread);
		}

		// then wait for them to be completed
		for (Thread pTaskThread : parallelTasks) {
			try {
				pTaskThread.join();
			} catch (Exception ignored) {
				System.err.println("Error from SimpleParallelTaskExecuter : "
						+ ignored.toString());
				ignored.printStackTrace();
			}
		}
	}

	/**
	 * Return the task list iterator
	 * 
	 * @return the task list
	 */
	public Iterator<SimpleParallelTask> getTaskList() {
		return pTaskList.iterator();
	}
}
