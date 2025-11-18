package transferlearn;

import allconstants.Constants;
import allconstants.NceConstants;
import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import peptideptmformatting.PeptideFormatter;
import utils.Print;
import utils.ProgressReporter;

import javax.net.ssl.HttpsURLConnection;
import java.io.*;
import java.lang.reflect.Type;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.*;
import java.util.HashMap;
import java.util.Map;

import static utils.Print.printInfo;

public class Helpers {
    static HttpURLConnection setUpConnection(String serverString, URL serverURL) throws IOException {
        //set up connection
        HttpURLConnection connection;
        if (serverString.startsWith("http:")) {
            connection = (HttpURLConnection) serverURL.openConnection();
        } else { //https:
            connection = (HttpsURLConnection) serverURL.openConnection();
        }
        return connection;
    }

    static HashMap<String, Object> readJsonResponse(InputStream stream) throws IOException {
        try (InputStreamReader reader = new InputStreamReader(stream, StandardCharsets.UTF_8)) {
            Gson gson = new Gson();

            // Use Object for values so nested objects are deserialized as Maps
            Type type = new TypeToken<HashMap<String, Object>>() {}.getType();

            return gson.fromJson(reader, type);
        }
    }

    static void convertCsvToParquet(String csvFilePath, String parquetFilePath, boolean deleteCsv) {
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {
            Statement stmt = conn.createStatement();

            // More control over CSV reading and Parquet writing
            String sql = String.format(
                    "COPY (SELECT * FROM read_csv('%s', " +
                            "  header=true, " +
                            "  auto_detect=true, " +
                            "  sample_size=1000" +
                            ")) TO '%s' (FORMAT PARQUET, " +
                            "  COMPRESSION 'SNAPPY', " +
                            "  ROW_GROUP_SIZE 100000" +
                            ")",
                    csvFilePath, parquetFilePath
            );

            stmt.execute(sql);
            Print.printInfo(csvFilePath + " successfully converted to parquet");
            if (deleteCsv) {
                Print.printInfo("Now deleting " + csvFilePath);
                File csvFile = new File(csvFilePath);
                csvFile.delete();
            }

        } catch (SQLException e) {
            Print.printError("Error converting " + csvFilePath + " to parquet: " + e.getMessage());
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    static void convertParquetToLibraryTsv(String parquet, String tsv, HashMap<String, String> protToGene)
            throws SQLException, IOException {
        Print.printInfo("Loading fasta to sql");
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {

            // Create temporary table with mapping
            stmt.execute("CREATE TEMP TABLE mapping (key VARCHAR, value VARCHAR)");

            // Batch insert mapping data (much faster than individual inserts)
            conn.setAutoCommit(false);
            try (PreparedStatement pstmt = conn.prepareStatement(
                    "INSERT INTO mapping VALUES (?, ?)")) {
                int count = 0;
                for (Map.Entry<String, String> entry : protToGene.entrySet()) {
                    pstmt.setString(1, entry.getKey());
                    pstmt.setString(2, entry.getValue());
                    pstmt.addBatch();

                    if (++count % 10000 == 0) {
                        pstmt.executeBatch();
                    }
                }
                pstmt.executeBatch();
            }
            conn.commit();
            conn.setAutoCommit(true);

            // Perform join and export in one SQL operation
            Print.printInfo("Converting parquet to library.tsv format");
            String query = String.format(
                "COPY (SELECT p.*, COALESCE(m.value, p.ProteinId) AS GeneName " +
                    "FROM (SELECT *, row_number() OVER () as row_num FROM read_parquet('%s')) p " +
                    "LEFT JOIN mapping m ON p.ProteinId = m.key " +
                    "ORDER BY p.row_num) " +
                    "TO '%s' (FORMAT CSV, DELIMITER '\t', HEADER);",
                parquet.replace("'", "''"),
                tsv.replace("'", "''")
            );
            stmt.execute(query);
        }
    }

    static void convertPeptideListToApdInput(String inputParquet, File outputParquet) {
        printInfo("Converting peptide list to AlphaPeptDeep parquet");

        Path tmpCsv = null;
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {

            // Create temp CSV file
            tmpCsv = Files.createTempFile("alphapeptdeep_", ".csv");

            // Prepare writer to CSV
            try (BufferedWriter writer = Files.newBufferedWriter(tmpCsv)) {
                writer.write("sequence,mods,mod_sites,charge,nce,instrument,modified,proteins,is_decoy\n");

                // Get number of lines for progress reporting
                ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM read_parquet('" + inputParquet + "')");
                rs.next();
                int lines = rs.getInt(1);
                ProgressReporter pr = new ProgressReporter(lines);
                Print.printInfo(lines + " entries to convert");

                // Query input parquet
                rs = stmt.executeQuery(
                        "SELECT peptide, proteins, is_decoy FROM read_parquet('" + inputParquet + "')");

                while (rs.next()) {
                    String peptide = rs.getString("peptide");
                    String proteins = rs.getString("proteins");
                    boolean isDecoy = rs.getBoolean("is_decoy");

                    // extract charge
                    String charge = "";
                    while (!peptide.isEmpty() && Character.isDigit(peptide.charAt(peptide.length() - 1))) {
                        charge = peptide.charAt(peptide.length() - 1) + charge;
                        peptide = peptide.substring(0, peptide.length() - 1);
                    }

                    PeptideFormatter pf = new PeptideFormatter(peptide, charge, "apdpred");

                    writer.write(String.join(",",
                            pf.getStripped(),
                            pf.getAlphapeptdeepMods(),
                            pf.getModPositions(),
                            charge,
                            String.valueOf(NceConstants.getNCE()),
                            Constants.instrument,
                            pf.getLibrarytsv(),
                            proteins,
                            String.valueOf(isDecoy)
                    ));
                    writer.write("\n");
                    pr.progress();
                }
            }

            stmt.execute("CREATE TABLE tmp AS SELECT * FROM read_csv_auto('" + tmpCsv + "')");
            stmt.execute("COPY tmp TO '" + outputParquet.getAbsolutePath() + "' (FORMAT PARQUET)");

        } catch (Exception e) {
            Print.printError("Error: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        } finally {
            if (tmpCsv != null) {
                try { Files.deleteIfExists(tmpCsv); } catch (IOException ignored) {}
            }
        }
    }

    static void convertPeptideListToCsv(String peptideList, File csvFile) {
        printInfo("Converting peptide list to AlphaPeptDeep format");
        try (BufferedReader reader = new BufferedReader(new FileReader(peptideList));
             BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile))
        ) {
            long lines = Files.lines(Paths.get(peptideList)).count() - 1;
            ProgressReporter pr = new ProgressReporter((int) lines);

            //skip old header, write new header
            String line = reader.readLine(); //peptide,proteins,is_decoy
            writer.write("sequence,mods,mod_sites,charge,nce,instrument,modified,proteins,is_decoy\n");
            while ((line = reader.readLine()) != null) {
                String[] lineSplit = line.split(",");

                String charge = "";
                while (Character.isDigit(lineSplit[0].charAt(lineSplit[0].length() - 1))) {
                    charge = lineSplit[0].charAt(lineSplit[0].length() - 1) + charge;
                    lineSplit[0] = lineSplit[0].substring(0, lineSplit[0].length() - 1);
                }

                PeptideFormatter pf = new PeptideFormatter(lineSplit[0], charge, "apdpred");
                if (charge.isEmpty()) {
                    writer.write(pf.getStripped() + "," + pf.getAlphapeptdeepMods() + "," +
                            pf.getModPositions() + "," + "," + NceConstants.getNCE() + "," +
                            Constants.instrument + "," + pf.getLibrarytsv() + "," +
                            lineSplit[1] + "," + lineSplit[2] + "\n");
                } else {
                    writer.write(pf.getStripped() + "," + pf.getAlphapeptdeepMods() + "," +
                            pf.getModPositions() + "," + charge + "," + NceConstants.getNCE() + "," +
                            Constants.instrument + "," + pf.getLibrarytsv() + "," +
                            lineSplit[1] + "," + lineSplit[2] + "\n");
                }
                pr.progress();
            }
        } catch (IOException e) {
            Print.printError("Error reading and writing input file for AlphaPeptDeep prediction: "
                    + e.getMessage());
            System.exit(1);
        }
    }

    static class EndJob extends Thread {
        private final String cancelUrlPath;

        public EndJob(String cancelUrlPath) {
            this.cancelUrlPath = cancelUrlPath;
        }
        public void run() {
            URL cancelUrl;
            try {
                cancelUrl = new URL(cancelUrlPath);
                HttpURLConnection connection = setUpConnection(cancelUrlPath, cancelUrl);
                connection.connect();
                InputStream responseStream = connection.getInputStream();
                BufferedReader in = new BufferedReader(new InputStreamReader(responseStream));
                String line;
                while ((line = in.readLine()) != null) {
                    System.out.println(line);
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    //use this method to deal with parquet column during prediction to get a map with no semicolon parsing needed in future
    public static HashMap<String, String> mapProteinsListToGenes(String parquetPath, HashMap<String, String> protToGene) {
        Print.printInfo("Mapping peptide list's proteins to genes");
        HashMap<String, String> finalMap = new HashMap<>();

        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {

            String query = String.format("SELECT proteins FROM read_parquet('%s')",
                    parquetPath.replace("'", "''"));

            try (ResultSet rs = stmt.executeQuery(query)) {
                while (rs.next()) {
                    String proteins = rs.getString("proteins");
                    if (proteins == null || proteins.isEmpty()) {
                        continue;
                    }

                    String gene = transformProteins(proteins, protToGene);
                    finalMap.put(proteins, gene);
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException(e);
        }

        Print.printInfo("Mapped " + finalMap.size() + " protein labels to gene labels");
        return finalMap;
    }

    private static String transformProteins(String proteins, HashMap<String, String> protToGene) {
        // Case 1: Single key exists in map
        if (protToGene.containsKey(proteins)) {
            String gene = protToGene.get(proteins);
            if (gene == null || gene.isEmpty()) {
                return extractFromPipeSplit(proteins);
            }
            return gene;
        }

        // Case 2: Semicolon-separated keys
        if (proteins.contains(";")) {
            String[] keys = proteins.split(";");
            StringBuilder result = new StringBuilder();

            for (int i = 0; i < keys.length; i++) {
                if (i > 0) result.append(";");
                String key = keys[i].trim();

                if (protToGene.containsKey(key)) {
                    result.append(protToGene.get(key));
                } else {
                    // Not in map, try pipe split
                    result.append(extractFromPipeSplit(key));
                }
            }
            return result.toString();
        }

        // Case 3: Not a key, split by "|" and return second element
        return extractFromPipeSplit(proteins);
    }

    private static String extractFromPipeSplit(String value) {
        if (value.contains("|")) {
            String[] parts = value.split("\\|");
            if (parts.length >= 2) {
                return parts[1];
            }
        }
        // If no pipe or not enough parts, return original value
        return value;
    }

}
