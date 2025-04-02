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

package utils;

import allconstants.Constants;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

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
    static HashSet<String> unispecModels = new HashSet<>(Arrays.asList("Lumos", "QE", "QEHFX", "ELITE", "VELOS", "NONE"));

    public static String getInstrument(ScanCollectionDefault scans) {
        if (Constants.instrument.isEmpty()) {
            try {
                scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
                String model = scans.getRunInfo().getDefaultInstrument().getModel();
                String analyzer = scans.getRunInfo().getDefaultInstrument().getAnalyzer();
                for (String k : LumosKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: Lumos");
                        Constants.instrument = "Lumos";
                        return "Lumos";
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
                        printInfo("Instrument detected: SciexTOF");
                        Constants.instrument = "SciexTOF";
                        return "SciexTOF";
                    }
                }
                for (String k : timsTOFKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: timsTOF");
                        Constants.instrument = "timsTOF";
                        return "timsTOF";
                    }
                }
                for (String k : ThermoTOFKeys) {
                    if (model.contains(k) || analyzer.contains(k)) {
                        printInfo("Instrument detected: ThermoTOF");
                        Constants.instrument = "ThermoTOF";
                        return "ThermoTOF";
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

    public static String mapInstrumentToModelSpecific(String model) {
        model = model.toLowerCase();
        switch(model){
            case "unispec":
                if (!unispecModels.contains(Constants.instrument)) {
                    return "NONE";
                }
                return Constants.instrument;
            default:
                return Constants.instrument;
        }
    }
}
