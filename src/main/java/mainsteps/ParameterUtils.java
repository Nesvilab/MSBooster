package mainsteps;

import allconstants.*;
import citations.KoinaYamlParser;
import features.spectra.MassCalculator;
import peptideptmformatting.PTMhandler;
import utils.Model;
import utils.NumericUtils;

import java.io.*;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import static peptideptmformatting.PTMhandler.tmtMasses;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class ParameterUtils {
    static void printHelpMessage() {
        System.out.println("Usage: java -jar MSBooster-1.2.57.jar [flags]");
        System.out.println("Usage: java -jar MSBooster-1.2.57.jar --paramsList {*.txt}");
        System.out.println("A tool for annotating PSM pin files with deep learning-based features. " +
                "This has been tested on MSFragger pin files and may generalize to pin files from any tool," +
                " although this has not been specifically tested. ");
        System.out.println("For more information on required and optional parameters, please refer to the " +
                "MSBooster Github at https://github.com/Nesvilab/MSBooster");
        System.out.println("Printing example paramsList to exampleParams.txt. It will be a bit verbose, so " +
                "please delete any parameters not relevant to your analysis.");
        printParams(".");
        System.exit(0);
    }

    public static void processCommandLineInputs(String[] args, HashMap<String, String> params) throws IOException {
        for (int i = 0; i < args.length; i++) {
            String key = args[i].substring(2); //remove --
            if (key.equals("help")) { //help message
                printHelpMessage();
            }
            i++;
            StringBuilder sb = new StringBuilder(args[i]);
            if (i + 1 >= args.length) {
                params.put(key, sb.toString());
            } else {
                while (!args[i + 1].startsWith("--")) {
                    i++;
                    sb.append(" ");
                    sb.append(args[i]);
                    if (i + 1 >= args.length) {
                        break;
                    }
                }
                params.put(key, sb.toString());
            }
            if (key.equals("paramsList")) {
                processParameterList(params);
            }
        }
    }

    static void processParameterList(HashMap<String, String> params) throws IOException {
        //params with nulls are left out
        String line;
        BufferedReader reader = new BufferedReader(new FileReader(params.get("paramsList")));
        while ((line = reader.readLine()) != null) {
            if (!line.contains("=")) { //maybe empty line or comment line with #
                continue;
            }
            String[] lineSplit = line.split("=", 2);

            //check if null here
            if (!lineSplit[1].trim().equals("null")) {
                params.put(lineSplit[0].trim(), lineSplit[1].trim());
            }
        }
        reader.close();
    }

    static void readFraggerParams(HashMap<String, String> params) throws IOException {
        String line;
        boolean ppmToDa = false; //set Da tolernace once all parameters read in. No guarantee which param is read first
        BufferedReader reader = new BufferedReader(new FileReader(params.get("fragger")));
        printInfo("Reading " + params.get("fragger"));
        while ((line = reader.readLine()) != null) {
            String[] lineSplit = line.split("#")[0].split("=");
            if (lineSplit.length != 2) {
                continue;
            }
            String key = lineSplit[0].trim();
            String val = lineSplit[1].trim();
            switch (key) {
                case "fragment_mass_tolerance":
                    params.put("ppmTolerance", val);
                    break;
                case "fragment_mass_units":
                    if (val.equals("0")) {
                        ppmToDa = true;
                    }
                    break;
                case "decoy_prefix":
                    params.put("decoyPrefix", ">" + val);
                    break;
                case "search_enzyme_cutafter":
                    params.put("cutAfter", val);
                    break;
                case "search_enzyme_butnotafter":
                    params.put("butNotAfter", val);
                    break;
                case "digest_min_length":
                    params.put("digestMinLength", val);
                    break;
                case "digest_max_length":
                    params.put("digestMaxLength", val);
                    break;
                case "digest_mass_range":
                    String[] vals = val.split(" ");
                    params.put("digestMinMass", vals[0]);
                    params.put("digestMaxMass", vals[1]);
                    break;
                case "database_name":
                    params.put("fasta", val);
                    break;
                case "precursor_charge":
                    vals = val.split(" ");
                    params.put("minPrecursorCharge",
                            String.valueOf(Math.min(Integer.parseInt(vals[0]),
                                    Constants.minPrecursorCharge)));
                    params.put("maxPrecursorCharge",
                            String.valueOf(Math.max(Integer.parseInt(vals[1]),
                                    Constants.maxPrecursorCharge)));
                    break;
                case "mass_offsets":
                    if (!val.isEmpty()) {
                        vals = val.split("/");
                        String final_vals = "";
                        for (String v : vals) {
                            float fv = Float.parseFloat(v);
                            if (fv != 0) {
                                if (!final_vals.isEmpty()) {
                                    final_vals += "&";
                                }
                                final_vals += String.format("%.4f", fv);
                            }
                        }
                        params.put("massOffsets", final_vals);
                    }
                    break;
                case "mass_offsets_detailed":
                    if (!val.isEmpty()) {
                        vals = val.split(";");
                        String final_vals = "";
                        for (String v : vals) {
                            v = v.split("\\(")[0];
                            float fv = Float.parseFloat(v);
                            if (fv != 0) {
                                if (!final_vals.isEmpty()) {
                                    final_vals += "&";
                                }
                                final_vals += v;
                            }
                        }
                        params.put("massOffsetsDetailed", final_vals);
                    }
                    break;
                case "mass_diff_to_variable_mod":
                    params.put("massDiffToVariableMod", val);
                    break;
                case "add_B_user_amino_acid":
                    if (! val.equals("0.0")) {
                        MassCalculator.AAmap.put('B', Float.valueOf(val));
                    }
                    break;
                case "add_J_user_amino_acid":
                    if (! val.equals("0.0")) {
                        MassCalculator.AAmap.put('J', Float.valueOf(val));
                    }
                    break;
                case "add_O_user_amino_acid":
                    if (! val.equals("0.0")) {
                        MassCalculator.AAmap.put('O', Float.valueOf(val));
                    }
                    break;
                case "add_U_user_amino_acid":
                    if (! val.equals("0.0")) {
                        MassCalculator.AAmap.put('U', Float.valueOf(val));
                    }
                    break;
                case "add_Z_user_amino_acid":
                    if (! val.equals("0.0")) {
                        MassCalculator.AAmap.put('Z', Float.valueOf(val));
                    }
                    break;
                case "add_X_user_amino_acid":
                    if (! val.equals("0.0")) {
                        MassCalculator.AAmap.put('X', Float.valueOf(val));
                    }
                    break;
                default:
                    if (key.startsWith("variable_mod_")) { //reading in TMT/iTraq values when variable
                        Double varmod = Double.valueOf(val.trim().split(" ")[0]);
                        for (double potentialTmtMass : tmtMasses) {
                            if (NumericUtils.massesCloseEnough(varmod, potentialTmtMass)) {
                                printInfo("TMT/iTRAQ mass detected in fragger.params as variable modification: " +
                                        potentialTmtMass);
                                PTMhandler.setTmtMass(potentialTmtMass);
                                ModelCollections.rtCollection = "isolabel";
                                ModelCollections.ms2Collection = "isolabel";
                                break;
                            }
                        }
                    } else if (key.startsWith("add_")) { //reading in TMT/iTraq values when fixed mode
                        for (double potentialTmtMass : tmtMasses) {
                            double fixedMod = Double.parseDouble(val);
                            if (NumericUtils.massesCloseEnough(fixedMod, potentialTmtMass)) {
                                printInfo("TMT/iTRAQ mass detected in fragger.params as fixed modification: " +
                                        potentialTmtMass);
                                PTMhandler.setTmtMass(potentialTmtMass);
                                ModelCollections.rtCollection = "isolabel";
                                ModelCollections.ms2Collection = "isolabel";
                                break;
                            }
                        }
                    }
            }
        }

        float tol = Float.parseFloat(params.get("ppmTolerance"));
        if (ppmToDa) { //read in from msfragger params. Low res tolerance used for all matching
            Constants.matchWithDaltons = true;
            Constants.matchWithDaltonsAux = true;
            Constants.matchWithDaltonsDefault = true;
            Constants.DaTolerance = tol;
            printInfo("Using Dalton tolerance of " + Constants.DaTolerance + " Da");
        } else {
            if (tol >= 100f) {
                params.put("lowResppmTolerance", String.valueOf(tol));
            } else {
                params.put("highResppmTolerance", String.valueOf(tol));
            }
        }
    }

    static void printParams(String directory) {
        try {
            Constants c = new Constants();
            BufferedWriter buffer = new BufferedWriter(new FileWriter(directory + File.separator + "exampleParams.txt"));

            Field[] f = Constants.class.getFields();
            for (Field field : f) {
                if ((field.getModifiers() & Modifier.FINAL) != Modifier.FINAL) {
                    if (!field.getName().equals("paramsList")) {
                        buffer.write(field.getName() + " = " + field.get(c) + "\n");
                    }
                }
            }
            buffer.close();
        } catch (Exception e) {
            printError("could not write final params");
            e.getStackTrace();
            System.exit(1);
        }
    }

    //for updating constants read from parameter file
    static void updateField(Field field, String value, ConstantsInterface c)
            throws NoSuchMethodException, InvocationTargetException, InstantiationException, IllegalAccessException {
        //get class of field
        Class<?> myClass = field.getType();

        //do not parse use[something] if null
        if (myClass.getTypeName().equals("java.lang.Boolean")) {
            if (value.equals("null")) {
                return;
            }
        }

        //parse to appropriate type
        field.set(c, myClass.getConstructor(String.class).newInstance(value));
    }

    public static void updateConstants(HashMap<String, String> params) throws Exception {
        HashSet<String> fields = new HashSet<>();
        for (Field f : Constants.class.getDeclaredFields()) {
            fields.add(f.getName());
        }
        Constants c = new Constants();

        HashSet<String> fieldsFragmentIon = new HashSet<>();
        for (Field f : FragmentIonConstants.class.getDeclaredFields()) {
            fieldsFragmentIon.add(f.getName());
        }
        FragmentIonConstants fic = new FragmentIonConstants();

        HashSet<String> fieldsNCE = new HashSet<>();
        for (Field f : NceConstants.class.getDeclaredFields()) {
            fieldsNCE.add(f.getName());
        }
        NceConstants nc = new NceConstants();

        HashSet<String> fieldsModelCollections = new HashSet<>();
        for (Field f : ModelCollections.class.getDeclaredFields()) {
            fieldsModelCollections.add(f.getName());
        }
        ModelCollections mc = new ModelCollections();

        for (Map.Entry<String, String> entry : params.entrySet()) {
            String key = entry.getKey();
            if (key.charAt(0) == '#') { //comment
                continue;
            }
            if (key.charAt(0) == '/') { //comment
                continue;
            }
            if (!fields.contains(key) && !fieldsFragmentIon.contains(key) &&
                    !fieldsNCE.contains(key) && !fieldsModelCollections.contains(key)) {
                throw new Exception(entry.getKey() + " is not a valid parameter");
            } else if (fields.contains(key)) {
                Field field = Constants.class.getField(key);
                updateField(field, entry.getValue(), c);
            } else if (fieldsFragmentIon.contains(key)) {
                Field field = FragmentIonConstants.class.getField(key);
                updateField(field, entry.getValue(), fic);
            } else if (fieldsNCE.contains(key)) {
                Field field = NceConstants.class.getField(key);
                updateField(field, entry.getValue(), nc);
            } else {
                Field field = ModelCollections.class.getField(key);
                updateField(field, entry.getValue(), mc);
            }
        }
        c.updateOutputDirectory();
        c.updateInputPaths(); //setting null paths
    }

    public static void finalizeParameterFile() throws IOException {
        finalizeParameterFile(new ArrayList<>());
    }

    static void finalizeParameterFile(ArrayList<Model> models) throws IOException {
        //modify parameter file with updated info
        // Read all lines from the file into a list
        List<String> lines = Files.readAllLines(Paths.get(Constants.paramsList));

        //remove specific lines that should not be repeated
        for (int l = lines.size() - 1; l > -1; l--) {
            if (lines.get(l).startsWith("#CITATIONS")) {
                lines.remove(l);
            }
        }

        AtomicBoolean spectrafound = new AtomicBoolean(false);
        AtomicBoolean rtfound = new AtomicBoolean(false);
        AtomicBoolean imfound = new AtomicBoolean(false);
        AtomicBoolean spectraPredFilePDVFound = new AtomicBoolean(false);

        // Stream the list, and if a line starts with the specified prefix,
        // replace everything after the prefix with newSuffix
        // Otherwise keep the line as is
        List<String> modifiedLines = lines.stream()
                .map(line -> {
                    if (line.trim().startsWith("spectraModel")) {
                        spectrafound.set(true);
                        return "spectraModel=" + Constants.spectraModel;
                    } else if (line.trim().startsWith("rtModel")) {
                        rtfound.set(true);
                        return "rtModel=" + Constants.rtModel;
                    } else if (line.trim().startsWith("imModel")) {
                        imfound.set(true);
                        return "imModel=" + Constants.imModel;
                    } else if (line.trim().startsWith("spectraPredFilePDV")) {
                        spectraPredFilePDVFound.set(true);
                        return "spectraPredFilePDV=" + Constants.spectraPredFile;
                    }
                    return line;
                })
                .collect(Collectors.toList());

        // If no line was found with the specified prefix, add the new line
        if (!spectrafound.get()) {
            modifiedLines.add("spectraModel=" + Constants.spectraModel);
        }
        if (!rtfound.get()) {
            modifiedLines.add("rtModel=" + Constants.rtModel);
        }
        if (!imfound.get()) {
            modifiedLines.add("imModel=" + Constants.imModel);
        }
        //FragPipe PDV needs to know which file has spectral predictions
        if (Constants.useSpectra && !spectraPredFilePDVFound.get()) {
            modifiedLines.add("spectraPredFilePDV=" + Constants.spectraPredFile);
        }

        // Write the modified lines back to the file
        Files.write(Paths.get(Constants.paramsList), modifiedLines);

        //write citations of models used
        if (!Constants.KoinaURL.isEmpty()) {
            KoinaYamlParser kyp = new KoinaYamlParser();
            kyp.writeCitations(models);
        }
    }
}
