package writers;

import readers.datareaders.PinReader;
import utils.Print;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;

public class LibraryTsvReducer {
    final String fullLib;
    public final String newLib;
    int sequenceIdx;
    int modifiedIdx;
    int chargeIdx;
    int maxIdx;
    final File[] pinFiles;
    public boolean needReduced = false;

    public LibraryTsvReducer(String fullLib, File[] pinFiles) {
        this.fullLib = fullLib;
        this.newLib = fullLib.replace(".tsv", "_reduced.tsv");
        this.pinFiles = pinFiles;
    }

    public LibraryTsvReducer(String fullLib, String[] pinStrings) {
        this.fullLib = fullLib;
        this.newLib = fullLib.replace(".tsv", "_reduced.tsv");
        File[] pinFiles = new File[pinStrings.length];
        for (int i = 0; i < pinStrings.length; i++) {
            pinFiles[i] = new File(pinStrings[i]);
        }
        this.pinFiles = pinFiles;
    }

    //NOTE: this is meant to be heuristic, so some peptides not in pin files may sneak through
    public void writeReducedLib() {
        try (BufferedReader reader = new BufferedReader(new FileReader(fullLib), 4 * 1024 * 1024);
             BufferedWriter writer = new BufferedWriter(new FileWriter(newLib), 4 * 1024 * 1024)) {

            //get indices of PeptideSequence and PrecursorCharge
            String line = reader.readLine();
            writer.write(line + "\n");
            String[] splits = line.split("\t");
            for (int i = 0; i < splits.length; i++) {
                switch (splits[i]) {
                    case "PeptideSequence":
                        sequenceIdx = i;
                        break;
                    case "ModifiedPeptideSequence":
                        modifiedIdx = i;
                        break;
                    case "PrecursorCharge":
                        chargeIdx = i;
                        break;
                }
            }
            maxIdx = Math.max(Math.max(sequenceIdx, chargeIdx), modifiedIdx);

            //given a list of pin files, collect stripped peptides and charges in a hashset
            HashSet<String> strippedChargeSet = new HashSet<>();
            for (File pinFile : pinFiles) {
                PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
                while (pinReader.next(true)) {
                    strippedChargeSet.add(pinReader.getPep().getStripped() +
                            pinReader.getPep().getCharge());
                }
            }

            StringBuilder builder = new StringBuilder(1024 * 1024);
            boolean printedMessage = false;
            line = reader.readLine();
            while (line != null) {
                String[] strings = getStrippedChargeAndModCharge(line);
                //group by precursor
                String modCharge = strings[1];

                //iterate through large library tsv, checking if stripped peptide and charge is in hashset
                String strippedCharge = strings[0];
                boolean captureLines = strippedChargeSet.contains(strippedCharge);

                if (captureLines) {
                    builder.append(line).append("\n");
                } else {
                    //full library contains a peptide not in reduced
                    needReduced = true;
                }
                while ((line = reader.readLine()) != null) {
                    //while precursor stays the same, continue onto next line
                    if (substringFromModifiedIdx(line).startsWith(modCharge)) {
                        if (captureLines) {
                            builder.append(line).append("\n");
                        }
                    } else {
                        break;
                    }
                }

                //write accepted entries to a new library tsv "reduced"
                if (captureLines) {
                    if (builder.length() > 900000) {  // ~900KB
                        //approximation for a file that only has precursors we need
                        //if 900KB comes and still no excess precursors
                        if (!needReduced) {
                            Print.printInfo("Reading full library");
                            reader.close();
                            writer.close();
                            Files.deleteIfExists(Paths.get(newLib));
                            return;
                        }
                        if (!printedMessage) {
                            Print.printInfo("Writing reduced library from " +
                                    new File(fullLib).getName() + " to " + new File(newLib).getName());
                            printedMessage = true;
                        }
                        writer.write(builder.toString());
                        builder.setLength(0);  // Clear for reuse
                    }
                }
            }

            // Flush remaining at end
            if (builder.length() > 0) {
                writer.write(builder.toString());
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private String[] getStrippedChargeAndModCharge(String s) {
        String stripped = "";
        String charge = "";
        String mod = "";

        int start = 0;
        for (int i = 0; i <= maxIdx; i++) {
            int end = s.indexOf("\t", start);
            if (i == sequenceIdx) {
                stripped = s.substring(start, end);
            } else if (i == chargeIdx) {
                charge = s.substring(start, end);
            } if (i == modifiedIdx) {
                mod = s.substring(start, end);
            }
            start = end + 1;
        }
        return new String[]{stripped + charge, mod + "\t" + charge};
    }

    //this assumes modified sequence column is followed immediately by charge
    private String substringFromModifiedIdx(String s) {
        int tabCount = 0;

        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '\t') {
                tabCount++;
                if (tabCount == modifiedIdx - 1) {
                    return s.substring(i + 1);
                }
            }
        }
        return "";
    }
}
