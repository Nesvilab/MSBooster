package Features;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class PredictionEntryHashMap extends ConcurrentHashMap<String, PredictionEntry> {
    public void filterTopFragments(ExecutorService executorService)
            throws ExecutionException, InterruptedException {
        String[] peptides = new String[this.size()];
        PredictionEntry[] predictions = new PredictionEntry[this.size()];
        int entryIdx = 0;
        for (Map.Entry<String, PredictionEntry> entry : this.entrySet()) {
            peptides[entryIdx] = entry.getKey();
            predictions[entryIdx] = entry.getValue();
            entryIdx++;
        }

        List<Future> futureList = new ArrayList<>(Constants.numThreads);
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (this.size() * (long) i) / Constants.numThreads;
            int end = (int) (this.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    PredictionEntry pe = predictions[j];
                    pe.filterFragments();
                    this.put(peptides[j], pe);
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }
}
