package citations;

import allconstants.Constants;
import allconstants.ModelCollections;
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

    private String getCitation(String model) {
        //naming changes needed for yaml
        model = model.replace("AlphaPept", "AlphaPeptDeep");

        String description = paths.get("/" + model + "/infer").get("post").get("description");
        description = description.split("### Citation")[1].
                replace("<br>", " ").
                replace("\r", "").
                replace("\n", "");
        return description;
    }

    public void writeCitations(ArrayList<Model> models) throws IOException {
        HashSet<String> citations = new HashSet<>();
        for (Model model : models) {
            String modelName = model.name;
            if (ModelCollections.KoinaModels.contains(modelName)) {
                String citation = getCitation(modelName);
                citations.add(citation);
            }
        }

        if (!citations.isEmpty()) {
            String starter = "If using models from Koina, please cite the following:";
            String koinaCitation = "Lautenbacher L, Yang KL, et al. " +
                    "Koina: Democratizing machine learning for proteomics research. " +
                    "bioRxiv. (2024) doi:10.1101/2024.06.01.596953";
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
