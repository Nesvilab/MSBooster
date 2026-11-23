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
import koinaclasses.KoinaModelSpec;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.lcmsrun.LCMSRunInfo;
import umich.ms.datatypes.scan.props.Instrument;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

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

    public static String getInstrument(ScanCollectionDefault scans) {
        if (Constants.instrument.isEmpty()) {
            try {
                Instrument instrument = defaultInstrument(scans);
                String model = instrument == null ? "" : Objects.toString(instrument.getModel(), "");
                String analyzer = instrument == null ? "" : Objects.toString(instrument.getAnalyzer(), "");
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
                printInfo(undetectedInstrumentMessage(model, analyzer));
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

    //the run's default instrument, or null when there is none to read (the scans of an MGF input are null)
    private static Instrument defaultInstrument(ScanCollectionDefault scans) throws FileParsingException {
        if (scans == null) {
            return null;
        }
        try {
            scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
            LCMSRunInfo runInfo = scans.getRunInfo();
            return runInfo == null ? null : runInfo.getDefaultInstrument();
        } catch (NullPointerException e) {
            //a reader failure on an odd run header; the QE default has always covered it
            return null;
        }
    }

    //names the model and analyzer that matched no keyword, so the log shows what the mzML said
    private static String undetectedInstrumentMessage(String model, String analyzer) {
        return "Could not detect instrument type from model [" + model + "], analyzer [" + analyzer + "]. " +
                "Setting to QE. " +
                "If a different instrument was used, specify using '--instrument' via the command line " +
                "or 'instrument=' in the param file.";
    }

    //the run's instrument as a Koina model is sent it (KoinaModelSpec.InstrumentMap)
    public static String mapInstrumentToModelSpecific(KoinaModelSpec spec, String charge) {
        return spec.instrumentMap().map(Constants.instrument, charge);
    }
}
