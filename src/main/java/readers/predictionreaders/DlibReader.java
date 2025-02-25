/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package readers.predictionreaders;

import allconstants.Constants;
import features.spectra.MassCalculator;
import mainsteps.PinMzmlMatcher;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.PinReader;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

import static utils.Print.printInfo;

//refer to https://bitbucket.org/searleb/encyclopedia/wiki/EncyclopeDIA%20File%20Formats
public class DlibReader implements LibraryPredictionMapper {
    ArrayList<String> filenames = new ArrayList<>();
    PredictionEntryHashMap allPreds = new PredictionEntryHashMap();

    public DlibReader(String dlib) throws SQLException, IOException, FileParsingException, ExecutionException, InterruptedException, URISyntaxException {
        File predsDirectory = new File(dlib);
        String[] predsFiles = predsDirectory.list();
        filenames = new ArrayList<String>();

        if (predsFiles == null) { //if user provided a file, not a directory
            filenames.add(dlib);
        } else { //user provides directory
            for (String predsFile : predsFiles) {
                if (predsFile.contains(".dlib")) {
                    filenames.add(dlib + File.separator + predsFile);
                }
            }
        }

        //read in original dlib
        printInfo("loading dlib");
        for (String dl : filenames) {
            Connection connection = DriverManager.getConnection("jdbc:sqlite:" + dl);
            Statement statement = connection.createStatement();
            try {
                String query = "select * from entries";
                ResultSet resultSet = statement.executeQuery(query);
                while (resultSet.next()) {
                    String charge = resultSet.getString("PrecursorCharge");
                    String name = resultSet.getString("PeptideSeq"); //all C is +57
                    name = name.replace("C",
                            "C[" + PTMhandler.carbamidomethylationMass + "]");
                    name = name + "|" + charge;
                    float RT = resultSet.getFloat("RTInSeconds");
                    float[] mzs = extractMassArray(resultSet.getBytes("MassArray"),
                            resultSet.getInt("MassEncodedLength"));
                    float[] intensities = extractIntensityArray(resultSet.getBytes("IntensityArray"),
                            resultSet.getInt("IntensityEncodedLength"));

                    //filter out zero intensity
                    int newSize = 0;
                    for (float f : intensities) {
                        if (f > 0) {
                            newSize += 1;
                        }
                    }
                    float[] finalMZs = new float[newSize];
                    float[] finalIntensities = new float[newSize];
                    int index = 0;
                    for (int i = 0; i < mzs.length; i++) {
                        if (intensities[i] > 0) {
                            finalMZs[index] = mzs[i];
                            finalIntensities[index] = intensities[i];
                            index += 1;
                        }
                    }

                    PredictionEntry newPred = new PredictionEntry(finalMZs, finalIntensities,
                            new int[0], new int[0], new String[0]);
                    newPred.setRT(RT);
                    allPreds.put(name, newPred);
                }
            } catch (SQLException | IOException | DataFormatException e) {
                e.printStackTrace();
                System.exit(1);
            } finally {
                statement.close();
                connection.close();
            }
        }

        //for pin in pin files
        printInfo("generating modified and decoy peptides");
        PinMzmlMatcher pmm = new PinMzmlMatcher(Constants.mzmlDirectory, Constants.pinPepXMLDirectory);
        for (File f : pmm.pinFiles) {
            PinReader pin = new PinReader(f.getCanonicalPath());

            //for row in pin
            while (pin.next(true)) {
                int td = pin.getTD();
                if (td == 0) { //decoy
                    PeptideFormatter pfDecoy = pin.getPep();

                    String startPeptide = pin.getRow()[pin.pepIdx];
                    String reversePeptide = "";
                    if (startPeptide.charAt(0) != '-') {
                        reversePeptide += startPeptide.charAt(0);
                    }
                    reversePeptide += startPeptide.substring(2, startPeptide.length() - 2);
                    boolean removedAA = false;
                    String targetPeptide = "";
                    for (int i = reversePeptide.length() - 1; i > -1; i--) {
                        char c = reversePeptide.charAt(i);

                        //need to be true for removedAA, letter, and capital
                        if ((! removedAA) && (Character.isAlphabetic(c)) &&
                                (Character.isUpperCase(c))) {
                            removedAA = true;
                        } else if ((removedAA) && (Character.isAlphabetic(c)) &&
                                (Character.isUpperCase(c))) {
                            targetPeptide += c;
                        }
                    }
                    targetPeptide = targetPeptide.replace("C",
                            "C[" + PTMhandler.carbamidomethylationMass + "]");
                    targetPeptide += "|" + pfDecoy.getCharge();

                    if (allPreds.containsKey(targetPeptide)) {
                        MassCalculator mcDecoy = new MassCalculator(pfDecoy.getBase(), pfDecoy.getCharge());
                        MassCalculator mcTarget = new MassCalculator(targetPeptide.split("\\|")[0],
                                pfDecoy.getCharge());
                        PredictionEntry pe = allPreds.get(targetPeptide);

                        String[] annotations = mcTarget.annotateMZs(pe.mzs, "default", false)[0];
                        float[] decoyMZs = new float[annotations.length];
                        for (int i = 0; i < annotations.length; i++) {
                            String anno = annotations[i];

                            String fragment = anno.split(";")[0];
                            String[] plusSplit = fragment.split("\\+");
                            float mz = mcDecoy.calcMass(Integer.parseInt(plusSplit[0].substring(1)),
                                    plusSplit[0].substring(0, 1), Integer.parseInt(plusSplit[1]), 0);
                            decoyMZs[i] = mz;
                        }
                        PredictionEntry newPred = new PredictionEntry(decoyMZs, pe.intensities,
                                new int[0], new int[0], new String[0]);
                        newPred.setRT(pe.RT);
                        allPreds.put(pfDecoy.getBaseCharge(), newPred);
                    }
                } else { //add modified peptides as needed
                    //if basecharge names are not the same, but stripped name is the same, add new entry
                    PeptideFormatter pf = pin.getPep();
                    if ((allPreds.containsKey(pf.getDlib())) && (!allPreds.containsKey(pf.getBaseCharge()))) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.getDlib());
                        MassCalculator mc = new MassCalculator(pf.getDlib().split("\\|")[0], pf.getCharge());
                        String[][] info = mc.annotateMZs(tmp.mzs, "default", false);
                        String[] annotations = info[0];

                        MassCalculator shiftedMC = new MassCalculator(pf.getBase(), pf.getCharge());
                        ArrayList<Float> newMZs = new ArrayList<>();
                        ArrayList<Float> newIntensities = new ArrayList<>();

                        for (int i = 0; i < annotations.length; i++) {
                            boolean assign = true;
                            //possbly multiple assignments
                            String[] colonSplit = annotations[i].split(";");
                            ArrayList<Float> finalmz = new ArrayList<>();
                            for (String anno : colonSplit) {
                                if (anno.equals("unknown")) {
                                    printInfo(pf.getDlib());
                                    printInfo(Arrays.toString(annotations));
                                    printInfo(Arrays.toString(tmp.mzs));
                                    printInfo(mc.fragmentIons.keySet().toString());
                                }

                                String[] plusSplit = anno.split("\\+");
                                float mz = shiftedMC.calcMass(Integer.parseInt(plusSplit[0].substring(1)),
                                        plusSplit[0].substring(0, 1), Integer.parseInt(plusSplit[1]), 0);

                                if (finalmz.size() == 0) {
                                    finalmz.add(mz);
                                } else {
                                    if (Math.abs(mz - finalmz.get(0)) < 0.1) {
                                        finalmz.add(mz);
                                    } else {
                                        assign = false;
                                        break;
                                    }
                                }
                            }
                            if (assign) {
                                newIntensities.add(tmp.intensities[i]);

                                float avg = 0f;
                                for (float fl : finalmz) {
                                    avg += fl;
                                }
                                newMZs.add(avg / finalmz.size());
                            }
                        }

                        float[] newMZsArray = new float[newMZs.size()];
                        float[] newIntensitiesArray = new float[newIntensities.size()];
                        for (int i = 0; i < newMZsArray.length; i++) {
                            newMZsArray[i] = newMZs.get(i);
                            newIntensitiesArray[i] = newIntensities.get(i);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry(newMZsArray, newIntensitiesArray,
                                new int[0], new int[0], new String[0]);
                        newPred.setRT(tmp.RT);
                        allPreds.put(pf.getBaseCharge(), newPred);
                    }
                }
            }
        }
    }

    private static float[] extractMassArray(byte[] compressedData, int uncompressedLength) throws IOException, DataFormatException {
        byte[] b = decompress(compressedData, uncompressedLength);
        double[] d=new double[b.length/8];
        ByteBuffer bb= ByteBuffer.wrap(b);
        bb.order(ByteOrder.BIG_ENDIAN);
        DoubleBuffer db=bb.asDoubleBuffer();
        db.get(d);

        float[] f = new float[d.length];
        for (int i = 0; i < f.length; i++) {
            f[i] = (float) d[i];
        }
        return f;
    }
    private static float[] extractIntensityArray(byte[] compressedData, int uncompressedLength) throws IOException, DataFormatException {
        byte[] b = decompress(compressedData, uncompressedLength);
        float[] d=new float[b.length/4];
        ByteBuffer bb= ByteBuffer.wrap(b);
        bb.order(ByteOrder.BIG_ENDIAN);
        FloatBuffer db=bb.asFloatBuffer();
        db.get(d);
        return d;
    }

    private static byte[] decompress(byte[] compressedData, int uncompressedLength) throws IOException, DataFormatException {
        Inflater decompresser=new Inflater();
        decompresser.setInput(compressedData);
        byte[] decompressedData=new byte[uncompressedLength];
        decompresser.inflate(decompressedData);
        decompresser.end();
        return decompressedData;
    }

    public PredictionEntryHashMap getPreds() throws IOException {
        return allPreds;
    }
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    public void clear() {
        allPreds.clear();
    }
}
