package citations;

import allconstants.Constants;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import org.yaml.snakeyaml.Yaml;
import utils.Model;
import utils.Print;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class KoinaYamlParser {
    private HashMap<String, HashMap<String, HashMap<String, String>>> paths;

    public KoinaYamlParser() {
        parse();
    }

    private void parse() {
        Yaml yaml = new Yaml();
        InputStream is = getClass().getClassLoader().getResourceAsStream("koina_citations.yml");
        HashMap<String, HashMap<String, HashMap<String, HashMap<String, String>>>> firstLayer
                = yaml.load(is);
        paths = firstLayer.get("paths");
    }

    //the citation under a model's path in koina_citations.yml (KoinaModelSpec.citationKey), blank if Koina gives
    //none, or null if the file has no such path
    String getCitation(String citationKey) {
        HashMap<String, HashMap<String, String>> path = paths.get("/" + citationKey + "/infer");
        if (path == null || path.get("post") == null || path.get("post").get("description") == null) {
            return null;
        }
        String[] description = path.get("post").get("description").split("### Citation");
        if (description.length < 2) {
            return null;
        }
        return description[1].
                replace("<br>", " ").
                replace("\r", "").
                replace("\n", "");
    }

    public void writeCitations(ArrayList<Model> models) throws IOException {
        HashSet<String> citations = new HashSet<>();
        for (Model model : models) {
            KoinaModelSpec spec = KoinaModelRegistry.find(model.name);
            if (spec != null) {
                String citation = getCitation(spec.citationKey());
                if (citation == null || citation.isBlank()) {
                    Print.printInfo("No citation found for Koina model " + spec.name() + " under "
                            + spec.citationKey() + " in koina_citations.yml");
                } else {
                    citations.add(citation);
                }
            }
        }

        if (!citations.isEmpty()) {
            String starter = "If using models from Koina, please cite the following:";
            //as koina_citations.yml (Koina's openapi.yml) asks Koina to be cited
            String koinaCitation = "Lautenbacher L, Yang KL, Kockmann T, et al. " +
                    "Koina: Democratizing machine learning for proteomics research. " +
                    "Nat Commun 16, 9933 (2025). doi:10.1038/s41467-025-64870-5";
            String comment = "#CITATIONS: ";

            //write to output
            Print.printInfo(starter);
            Print.printInfo(koinaCitation);
            for (String citation : citations) {
                Print.printInfo(citation);
            }

            //write to msbooster parameter file
            ArrayList<String> citationList = new ArrayList<>();
            citationList.add(comment + starter);
            citationList.add(comment + koinaCitation);
            for (String citation : citations) {
                citationList.add(comment + citation);
            }

            Files.write(Paths.get(Constants.paramsList), citationList, StandardOpenOption.APPEND);
        }
    }

    public static void main(String[] args) {
        KoinaYamlParser yamlParser = new KoinaYamlParser();
        for (String model : yamlParser.paths.keySet()) {
            model = model.split("/")[1];
            System.out.println(model);
            System.out.println(yamlParser.getCitation(model));
        }
    }
}
