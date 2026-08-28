package transferlearn;

import allconstants.Constants;
import allconstants.NceConstants;
import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import peptideptmformatting.PeptideFormatter;
import utils.Multithreader;
import utils.Print;
import utils.ProgressReporter;
import utils.ProteinLabels;

import javax.net.ssl.HttpsURLConnection;
import java.io.*;
import java.lang.reflect.Type;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.sql.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static utils.Print.printInfo;

public class Helpers {
    /** The column naming a peptide's leading protein, which every library carries. */
    private static final String PROTEIN_ID = "ProteinId";

    /** The column naming a peptide's whole protein group; only FragCast returns one. */
    private static final String ALL_MAPPED_PROTEINS = "AllMappedProteins";

    /**
     * That group reduced to the form the gene map is keyed in: {@code sp|P1|A_HUMAN;sp|P2|B_HUMAN}
     * becomes {@code P1;P2}.
     *
     * <p>This has to reproduce {@link #transformProteins} exactly, key for key, or the lookup misses
     * and the genes fall back to accessions. So it mirrors both of that method's rules: an accession
     * is the second pipe-delimited token when there is one and the label itself otherwise (which is
     * {@link #extractFromPipeSplit}), and a decoy keeps its prefix - MSFragger writes that prefix on
     * the {@code sp} token, so taking the accession alone would silently turn a decoy into its
     * target.
     *
     * <p>What it does not reproduce is the de-duplication {@link ProteinLabels#normalize} does, and it
     * does not have to: every writer normalises before a predictor ever sees the column, so a group
     * reaching a library cannot repeat a protein. A group that somehow did would miss this lookup and
     * fall back to its accessions, which is a duller answer rather than a wrong one.
     */
    private static String accessionsOfGroup() {
        final String prefix = Constants.decoyPrefix.replace("'", "''");
        return "array_to_string(list_transform(string_split(p." + ALL_MAPPED_PROTEINS + ", ';'), " +
                "x -> CASE WHEN len(string_split(x, '|')) >= 2 " +
                "THEN (CASE WHEN starts_with(x, '" + prefix + "') THEN '" + prefix + "' ELSE '' END) " +
                "|| string_split(x, '|')[2] " +
                "ELSE x END), ';')";
    }

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

    static void convertParquetToLibraryTsv(String parquet, String tsv, HashMap<String, String> protToGene,
                                           Connection conn) throws SQLException, IOException {
        try (Statement stmt = conn.createStatement()) {

            stmt.execute("CREATE TABLE IF NOT EXISTS mapping (key VARCHAR, value VARCHAR)");

            try (ResultSet check = stmt.executeQuery("SELECT COUNT(*) FROM mapping")) {
                check.next();
                if (check.getLong(1) == 0) {
                    Print.printInfo("Loading fasta to sql");
                    conn.setAutoCommit(false);
                    try (PreparedStatement pstmt = conn.prepareStatement("INSERT INTO mapping VALUES (?, ?)")) {
                        int count = 0;
                        for (Map.Entry<String, String> entry : protToGene.entrySet()) {
                            pstmt.setString(1, entry.getKey());
                            pstmt.setString(2, entry.getValue());
                            pstmt.addBatch();
                            if (++count % 10000 == 0) pstmt.executeBatch();
                        }
                        pstmt.executeBatch();
                    }
                    conn.commit();
                    conn.setAutoCommit(true);
                    stmt.execute("CREATE INDEX IF NOT EXISTS mapping_idx ON mapping(key)");
                }
            }

            Print.printInfo("Converting parquet to library.tsv format");

            // Fetch column names from parquet
            try (ResultSet rs = stmt.executeQuery(
                    String.format("SELECT * FROM read_parquet('%s') LIMIT 1", parquet.replace("'", "''")))) {

                ResultSetMetaData meta = rs.getMetaData();
                int colCount = meta.getColumnCount();

                //Only FragCast returns a protein group; the server has no such column, and then there
                //is no group to name genes for. Whether this library has one falls out of building
                //the select, so it is not looked for twice.
                final String groupKey = accessionsOfGroup();
                boolean hasGroup = false;

                //Each gene column follows the protein column it names, and is derived here rather than
                //kept from the predictor: a predictor reads its gene off the protein header, which
                //gives the entry mnemonic and not the gene. Both come from the FASTA instead.
                StringBuilder select = new StringBuilder();
                for (int i = 1; i <= colCount; i++) {
                    String col = meta.getColumnName(i);
                    select.append("p.").append(col);
                    if (PROTEIN_ID.equals(col)) {
                        select.append(", COALESCE(m.value, p.").append(PROTEIN_ID).append(") AS GeneName");
                    } else if (ALL_MAPPED_PROTEINS.equals(col)) {
                        hasGroup = true;
                        select.append(", COALESCE(g.value, ").append(groupKey).append(") AS AllMappedGenes");
                    }
                    if (i != colCount) select.append(", ");
                }

                // Write header if file is new/empty
                File tsvFile = new File(tsv);
                if (!tsvFile.exists() || tsvFile.length() == 0) {
                    try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(tsv),
                            StandardCharsets.UTF_8, StandardOpenOption.CREATE)) {
                        //Built from the same two rules the select is, so the two cannot fall out of
                        //step and shift every column of the library one place against its name.
                        StringBuilder header = new StringBuilder();
                        for (int i = 1; i <= colCount; i++) {
                            String col = meta.getColumnName(i);
                            if (i > 1) header.append('\t');
                            header.append(col);
                            if (PROTEIN_ID.equals(col)) {
                                header.append('\t').append("GeneName");
                            } else if (ALL_MAPPED_PROTEINS.equals(col)) {
                                header.append('\t').append("AllMappedGenes");
                            }
                        }
                        writer.write(header.toString());
                        writer.newLine();
                    }
                }

                // Write to temp file using native DuckDB COPY (no JDBC row iteration)
                Path tmpTsv = Files.createTempFile("apd_chunk_", ".tsv");
                try {
                    //The group's genes are looked up under the group's own key, which the map holds
                    //alongside the per-protein ones - the accessions in the order the predictor was
                    //given them, which is the order it wrote AllMappedProteins in.
                    String groupJoin = hasGroup ? " LEFT JOIN mapping g ON g.key = " + groupKey : "";
                    String query = String.format(
                            "COPY (" +
                                    "SELECT %s " +
                                    "FROM read_parquet('%s') p " +
                                    "LEFT JOIN mapping m ON p.ProteinId = m.key%s" +
                                    ") TO '%s' (FORMAT CSV, DELIMITER '\t', HEADER false);",
                            select,
                            parquet.replace("'", "''"),
                            groupJoin,
                            tmpTsv.toString().replace("\\", "/")
                    );
                    stmt.execute(query);

                    // Append temp file to TSV
                    try (OutputStream out = Files.newOutputStream(Paths.get(tsv),
                            StandardOpenOption.APPEND, StandardOpenOption.CREATE)) {
                        Files.copy(tmpTsv, out);
                    }
                } finally {
                    Files.deleteIfExists(tmpTsv);
                }
            }
        }
    }

    static File[] convertPeptideListToApdInput(String inputParquet, String outputParquetBasename, int maxSize){
        File[] inputFiles = null;
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            // Get number of lines for progress reporting and splitting into smaller files
            ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM read_parquet('" + inputParquet + "')");
            rs.next();
            int lines = rs.getInt(1);
            ProgressReporter pr = new ProgressReporter(lines);
            Print.printInfo(lines + " peptides to convert");

            //split into even groups of max size maxSize
            Multithreader mt = new Multithreader(lines, (int) Math.ceil((double) lines / (double) maxSize));
            inputFiles = new File[mt.indices.length - 1];

            // Query input parquet
            rs = stmt.executeQuery("SELECT peptide, proteins, is_decoy FROM read_parquet('" + inputParquet + "')");

            //send to mini method
            printInfo("Converting peptide list to AlphaPeptDeep parquet");
            for (int i = 0; i < mt.indices.length - 1; i++) {
                File inputFile = convertPeptideListToApdInput(rs, outputParquetBasename + i,
                        pr, mt.indices[i + 1] - mt.indices[i]);
                inputFiles[i] = inputFile;
            }
        } catch (SQLException e) {
            throw new RuntimeException(e);
        }
        return inputFiles;
    }

    static File convertPeptideListToApdInput(ResultSet rs, String outputParquetBasename,
                                             ProgressReporter pr, int numIter) {
        Path tmpCsv = null;
        String outputParquetPath = outputParquetBasename + ".parquet";
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {

            // Create temp CSV file
            tmpCsv = Files.createTempFile("alphapeptdeep_", ".csv");

            // Prepare writer to CSV
            try (BufferedWriter writer = Files.newBufferedWriter(tmpCsv)) {
                writer.write("sequence,mods,mod_sites,charge,nce,instrument,modified,proteins,is_decoy\n");

                int i = 0;
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
                            //normalised like every other protein column MSBooster writes: the server
                            //reads a ';' here as "this peptide is shared" (see ProteinLabels), and a
                            //null column would otherwise be joined in as the text "null"
                            ProteinLabels.normalize(proteins),
                            String.valueOf(isDecoy)
                    ));
                    writer.write("\n");
                    pr.progress();
                    i += 1;
                    if (i >= numIter) {
                        break;
                    }
                }
            }

            stmt.execute("CREATE TABLE tmp AS SELECT * FROM read_csv_auto('" + tmpCsv + "')");
            stmt.execute("COPY tmp TO '" + outputParquetPath + "' (FORMAT PARQUET)");
            return new File(outputParquetPath);
        } catch (Exception e) {
            Print.printError("Error: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        } finally {
            if (tmpCsv != null) {
                try { Files.deleteIfExists(tmpCsv); } catch (IOException ignored) {}
            }
        }
        return new File(outputParquetPath);
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
        public boolean ended = false;

        public EndJob(String cancelUrlPath) {
            this.cancelUrlPath = cancelUrlPath;
        }
        public void run() {
            if (!ended) {
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
    }

    /**
     * Map each peptide list's protein labels to gene labels, so no semicolon parsing is needed later.
     *
     * <p>The list is read as Parquet or as delimited text, by extension. The server path always has
     * Parquet; the local one writes its peptides as TSV for FragCast to read, and that file carries
     * the same {@code proteins} column, so it is the same question asked of a different container.
     */
    public static HashMap<String, String> mapProteinsListToGenes(String peptideListPath, HashMap<String, String> protToGene) {
        Print.printInfo("Mapping peptide list's proteins to genes");
        HashMap<String, String> finalMap = new HashMap<>();

        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {

            String lower = peptideListPath.toLowerCase();
            //The delimiter is named rather than sniffed. A peptide list's protein column is full of
            //pipes, and on a short list a sniffer has more of those to go on than it has delimiters -
            //it then reads 'sp|P1|A_HUMAN' as three columns, and the run dies on "Referenced column
            //proteins not found". Only small jobs are affected, which is the worst size to fail at:
            //nothing about a 5-peptide list suggests it should behave differently from a 50000-peptide
            //one. The extension already says which delimiter it is.
            String reader = (lower.endsWith(".parquet") || lower.endsWith(".pq"))
                    ? "read_parquet('%s')"
                    : "read_csv('%s', header=true, auto_detect=true, delim='"
                        + (lower.endsWith(".csv") ? "," : "\\t") + "')";
            String query = String.format("SELECT proteins FROM " + reader,
                    peptideListPath.replace("'", "''"));

            try (ResultSet rs = stmt.executeQuery(query)) {
                while (rs.next()) {
                    //normalised for the same reason the writers normalise: this map is keyed by the
                    //protein label the predictor was given, so reading the column by a different rule
                    //than the writers would key it on something the library never contains
                    String proteins = ProteinLabels.normalize(rs.getString("proteins"));
                    if (proteins.isEmpty()) {
                        continue;
                    }

                    String[] shorthandProtein_gene = transformProteins(proteins, protToGene);
                    finalMap.put(shorthandProtein_gene[0], shorthandProtein_gene[1]);
                    addEachProtein(finalMap, shorthandProtein_gene[0], shorthandProtein_gene[1]);
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException(e);
        }

        Print.printInfo("Mapped " + finalMap.size() + " protein labels to gene labels");
        return finalMap;
    }

    /**
     * Key a shared peptide's group by each of its proteins as well as by the group as a whole.
     *
     * <p>The group entry alone answers only the question the AlphaPeptDeep server asks. Both other
     * lookups are per protein: FragCast names a library's {@code ProteinId} after the FIRST protein
     * of the group, so the librarytsv join sees {@code P1} where the map held {@code P1;P2}, and
     * {@code ParquetToSpecLib} splits {@code ProteinId} and looks its parts up one at a time
     * whichever backend wrote it. Both miss, and a gene name silently degrades to the accession it
     * was supposed to replace.
     *
     * <p>The two strings are built in lockstep by {@link #transformProteins}, so their parts pair up
     * by position; a length that disagrees means something was not built that way, and is left alone
     * rather than paired up wrongly. A protein the list also carries on its own keeps that entry -
     * {@code putIfAbsent} never displaces one, and the standalone row's own {@code put} overwrites
     * whenever it is read - so nothing a group implies can outrank what a peptide states directly.
     */
    private static void addEachProtein(HashMap<String, String> finalMap, String shorthand, String genes) {
        if (!shorthand.contains(";")) {
            return;
        }
        String[] ids = shorthand.split(";");
        String[] geneNames = genes.split(";");
        if (ids.length != geneNames.length) {
            return;
        }
        for (int i = 0; i < ids.length; i++) {
            finalMap.putIfAbsent(ids[i], geneNames[i]);
        }
    }

    /**
     * A peptide's protein labels as the pair {shorthand accessions, gene names}, ';'-joined.
     *
     * <p>One protein and several are the same code deliberately. They used to be two branches, and
     * the branches drifted: only the several-protein one carried the decoy prefix over, so a peptide
     * sitting in a single decoy was keyed - and named - as though it were the target that decoy was
     * reversed from. Both predictors label a decoy by prefixing the accession they read out of its
     * header, so that key is one nothing in a library ever matches. Splitting on ';' yields one
     * element when there is no ';', which makes the two cases the same question.
     */
    private static String[] transformProteins(String proteins, HashMap<String, String> protToGene) {
        StringBuilder shorthandProteins = new StringBuilder();
        StringBuilder result = new StringBuilder();
        String[] labels = proteins.split(";");
        for (int i = 0; i < labels.length; i++) {
            if (i > 0) {
                shorthandProteins.append(";");
                result.append(";");
            }
            String label = labels[i].trim();
            String key = extractFromPipeSplit(label);
            //The prefix sits on the database token (rev_sp|P1|A_HUMAN), so an accession read out of
            //the header has lost it and has to have it put back. Only then: a label with nothing to
            //extract is returned whole, prefix included, and prefixing that again gives rev_rev_P1.
            String marker = !key.equals(label) && !Constants.decoyPrefix.isEmpty()
                    && label.startsWith(Constants.decoyPrefix)
                    ? Constants.decoyPrefix
                    : "";
            String gene = protToGene.get(key);
            shorthandProteins.append(marker).append(key);
            //A protein the FASTA does not name, or names with no gene, answers with its own accession
            //- the label is still right, it just says no more than the accession did.
            result.append(marker).append(gene == null || gene.isEmpty() ? key : gene);
        }
        return new String[]{shorthandProteins.toString(), result.toString()};
    }

    public static String extractFromPipeSplit(String value) {
        if (value.contains("|")) {
            String[] parts = value.split("\\|");
            if (parts.length >= 2) {
                return parts[1];
            }
        }
        // If no pipe or not enough parts, return original value
        return value;
    }

    public static void customModsStringToTsv(String mods) {
        //hardcoded to write to path
        try {
            Path tempFile = Files.createTempFile("custom_mods_", ".tsv");
            tempFile.toFile().deleteOnExit();
            Constants.additionalModsFile = tempFile.toFile().getAbsolutePath();
            BufferedWriter writer = new BufferedWriter(new FileWriter(Constants.additionalModsFile));
            writer.write("mod_name\tcomposition\tmodloss_composition\tmod_mass\n");

            String[] individualMods = mods.split(";");
            for (String mod : individualMods) {
                String line = mod.replace(",", "\t");
                writer.write(line + "\n");
            }
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
