package utils;

public class Multithreader {
    public int[] indices;
    //create make an array of equally sized entry subsets
    public Multithreader(int entries, int divisions) {
        indices = new int[divisions + 1];
        int subsetSize = entries / divisions;
        int remainder = entries % divisions;
        indices[0] = 0;
        for (int i = 1; i < divisions; i++) {
            indices[i] = indices[i - 1] + subsetSize + (i - 1 < remainder ? 1 : 0);
        }
        indices[divisions] = entries;
    }
}
