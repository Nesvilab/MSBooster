/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package utils;

import allconstants.Constants;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static utils.Print.printError;
import static utils.Print.printInfo;

public class InstrumentUtils {
    //TODO: support for astral model?

    //list of tags in mzml metadata that help ID instrument
    static HashSet<String> LumosKeys = new HashSet<>(Arrays.asList("LTQ", "Lumos", "Fusion", "Elite", "Velos", "Eclipse", "Tribrid"));
    static HashSet<String> QEKeys = new HashSet<>(Arrays.asList("QE", "Exactive", "Exploris"));
    static HashSet<String> SciexTOFKeys = new HashSet<>(Arrays.asList("Sciex", "TripleTOF"));
    static HashSet<String> timsTOFKeys = new HashSet<>(List.of("flight", "timsTOF"));
    static HashSet<String> ThermoTOFKeys = new HashSet<>(List.of("Astral"));

    //supported instruments by each model
    static CaseInsensitiveHashSet unispecInstruments = new CaseInsensitiveHashSet(Arrays.asList(
            "LUMOS", "QE", "QEHFX", "ELITE", "VELOS", "NONE"));
    static CaseInsensitiveHashSet alphapeptdeepInstruments = new CaseInsensitiveHashSet(Arrays.asList(
            "QE", "LUMOS", "TIMSTOF", "SCIEXTOF")); //"THERMOTOF"
//    static HashSet<String> allInstruments = makeAllInstruments();
//    private static HashSet<String> makeAllInstruments() {
//        HashSet<String> allInstruments = new HashSet<>();
//        allInstruments.addAll(unispecInstruments);
//        allInstruments.addAll(alphapeptdeepInstruments);
//        return allInstruments;
//    }

    public static String getInstrument(ScanCollectionDefault scans) {
        if (Constants.instrument.isEmpty()) {
            try {
                scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
                String model = scans.getRunInfo().getDefaultInstrument().getModel();
                String analyzer = scans.getRunInfo().getDefaultInstrument().getAnalyzer();
                for (String k : LumosKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: LUMOS");
                        Constants.instrument = "LUMOS";
                        return "LUMOS";
                    }
                }
                for (String k : QEKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: QE");
                        Constants.instrument = "QE";
                        return "QE";
                    }
                }
                for (String k : SciexTOFKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: SCIEXTOF");
                        Constants.instrument = "SCIEXTOF";
                        return "SCIEXTOF";
                    }
                }
                for (String k : timsTOFKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: TIMSTOF");
                        Constants.instrument = "TIMSTOF";
                        return "TIMSTOF";
                    }
                }
                for (String k : ThermoTOFKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: Astral");
                        Constants.instrument = "THERMOTOF";
                        return "THERMOTOF";
                    }
                }
                printInfo("Could not detect instrument type. Setting to QE. " +
                        "If a different instrument was used, specify using '--instrument' via the command line " +
                        "or 'instrument=' in the param file.");
                Constants.instrument = "QE";
                return "QE"; //default if nothing found
            } catch (NullPointerException e) {
                printInfo("Could not detect instrument type. Setting to QE. " +
                        "If a different instrument was used, specify using '--instrument' via the command line " +
                        "or 'instrument=' in the param file.");
                Constants.instrument = "QE";
                return "QE"; //default if nothing found
            } catch (FileParsingException e) {
                e.printStackTrace();
                System.exit(1);
                return "";
            }
        } else {
            return Constants.instrument;
        }
    }

    public static String mapInstrumentToModelSpecific(String model, String charge) {
        model = model.toLowerCase();
        String instrument = Constants.instrument;

        if (model.contains("unispec")) {
            if (!unispecInstruments.contains(instrument)) {
                instrument = "NONE";
            } else if ((instrument.equals("QE") || instrument.equals("QEHFX")) &&
                    charge.equals("1")) {
                instrument = "NONE";
            }
        } else if (model.contains("alphapept")) {
            if (instrument.equals("QEHFX")) {
                instrument = "QE";
            } else if (!alphapeptdeepInstruments.contains(instrument)) {
                instrument = "QE";
            }
        }

        return instrument;
    }
}
