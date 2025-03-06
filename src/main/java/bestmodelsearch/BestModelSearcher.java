package bestmodelsearch;

import utils.StatMethods;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import static utils.Print.printInfo;

public class BestModelSearcher {
    HashMap<String, float[]> scoreMap = new HashMap<>();
    LinkedHashMap<String, Float> summaryMap = new LinkedHashMap<>();
    public BestModelSearcher() {}

    private void makeSameLength(HashMap<String, float[]> scores, boolean biggerIsBetter) {
        int numEntries = Integer.MAX_VALUE;
        for (float[] diffs : scores.values()) {
            if (diffs.length < numEntries) {
                numEntries = diffs.length;
            }
        }

        for (Map.Entry<String, float[]> entry : scores.entrySet()) {
            float[] array = entry.getValue();
            Arrays.sort(array);
            float[] newArray;
            if (biggerIsBetter) {
                newArray = Arrays.copyOfRange(array, array.length - numEntries, array.length);
            } else {
                newArray = Arrays.copyOfRange(array, 0, numEntries);
            }
            scoreMap.put(entry.getKey(), newArray);
        }
    }

    private void makeSummaryHashMap(HashMap<String, Float> map) {
        map.entrySet().stream().sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
                .forEachOrdered(entry -> summaryMap.put(entry.getKey(), entry.getValue()));
    }

    private String selectModel(boolean biggerIsBetter) {
        float bestScore;
        if (biggerIsBetter) {
            bestScore = Float.MIN_VALUE;
        } else {
            bestScore = Float.MAX_VALUE;
        }

        String bestModel = "";
        for (Map.Entry<String, Float> entry : summaryMap.entrySet()) {
            if (biggerIsBetter) {
                if (entry.getValue() > bestScore) {
                    bestModel = entry.getKey();
                    bestScore = entry.getValue();
                }
            } else {
                if (entry.getValue() < bestScore) {
                    bestModel = entry.getKey();
                    bestScore = entry.getValue();
                }
            }
        }

        return bestModel;
    }

    //top: consider largest values (e.g. largest deltaRTs or highest spectral similarities)
    //bottom: opposite
    public String consensus(HashMap<String, float[]> scores, int topN, boolean biggerIsBetter, boolean top) {
        HashMap<String, Float> votes = new HashMap<>();
        for (String model : scores.keySet()) {
            votes.put(model, 0f);
        }

        makeSameLength(scores, biggerIsBetter);

        for (int i = 0; i < topN; i++) {
            String bestModel = "";
            float bestDiff;
            if (biggerIsBetter) {
                bestDiff = Float.MIN_VALUE;
            } else {
                bestDiff = Float.MAX_VALUE;
            }

            for (String model : scoreMap.keySet()) {
                float[] diffs = scoreMap.get(model);
                float delta;
                if (top) {
                    delta = diffs[diffs.length - 1 - i];
                } else {
                    delta = diffs[i];
                }
                if (biggerIsBetter) {
                    if (delta > bestDiff) {
                        bestDiff = delta;
                        bestModel = model;
                    }
                } else {
                    if (delta < bestDiff) {
                        bestDiff = delta;
                        bestModel = model;
                    }
                }
            }

            votes.put(bestModel, votes.get(bestModel) + 1f);
        }
        makeSummaryHashMap(votes);
        for (Map.Entry<String, Float> entry : summaryMap.entrySet()) {
            printInfo("Votes for " + entry.getKey() + ": " + String.format("%.0f", entry.getValue()));
        }
        return selectModel(true); //true because we want the most votes
    }

    public String medianMethod(HashMap<String, float[]> scores, boolean biggerIsBetter) {
        HashMap<String, Float> medians = new HashMap<>();
        makeSameLength(scores, biggerIsBetter);

        for (Map.Entry<String, float[]> entry : scoreMap.entrySet()) {
            double median = StatMethods.median(entry.getValue());
            medians.put(entry.getKey(), (float) median);
        }
        makeSummaryHashMap(medians);
        for (Map.Entry<String, Float> entry : summaryMap.entrySet()) {
            printInfo("Median similarity for " + entry.getKey() + " is " + String.format("%.4f", entry.getValue()));
        }
        return selectModel(biggerIsBetter);
    }

    public String RMSE(HashMap<String, float[]> scores, boolean biggerIsBetter) {
        HashMap<String, Float> rmses = new HashMap<>();
        makeSameLength(scores, biggerIsBetter);

        for (Map.Entry<String, float[]> entry : scoreMap.entrySet()) {
            float rmse = (float) Math.sqrt(StatMethods.meanSquaredError(entry.getValue()));
            rmses.put(entry.getKey(), rmse);
        }
        makeSummaryHashMap(rmses);
        for (Map.Entry<String, Float> entry : summaryMap.entrySet()) {
            printInfo(entry.getKey() + " has root mean squared error of " + String.format("%.4f", entry.getValue()));
        }
        return selectModel(biggerIsBetter);
    }

    public void writeScores(String outputFile) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
        StringBuilder sb = new StringBuilder();
        ArrayList<String> models = new ArrayList<>();
        for (String model : scoreMap.keySet()) {
            sb.append(model).append("\t");
            models.add(model);
        }
        sb.deleteCharAt(sb.length() - 1);
        bw.write(sb + "\n");
        
        int rowsToWrite = 0;
        for (Map.Entry<String, float[]> entry : scoreMap.entrySet()) {
            rowsToWrite = entry.getValue().length;
            break;
        }
        
        for (int i = 0; i < rowsToWrite; i++) {
            for (String model : models) {
                bw.write(scoreMap.get(model)[i] + "\t");
            }
            bw.write("\n");
        }
        bw.close();
    }
}
