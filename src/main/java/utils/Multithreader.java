package utils;

public class Multithreader {
    public int[] indices;
    //create make an array of equally sized entry subsets
    public Multithreader(int entries, int divisions) {
        indices = new int[divisions + 1];
        int subsetSize = (int) (entries / (float) divisions);
        for (int i = 0; i < divisions; i++) {
            indices[i] = subsetSize * i;
        }
        indices[divisions] = entries;
    }
}
