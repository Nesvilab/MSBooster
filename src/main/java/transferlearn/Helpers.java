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
import java.nio.file.Paths;
import java.sql.*;
import java.util.HashMap;

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

    static void convertParquetToCsv(String parquet, String tsv) throws SQLException, IOException {
        Print.printInfo("Converting parquet to CSV");
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             BufferedWriter writer = new BufferedWriter(new FileWriter(tsv))) {

            //get number of rows
            ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM read_parquet('" + parquet + "')");
            rs.next();
            int rowCount = rs.getInt(1);
            ProgressReporter pr = new ProgressReporter(rowCount);

            // Register and query Parquet
            String query = "SELECT * FROM read_parquet('" + parquet + "')";
            rs = stmt.executeQuery(query);

            // Write header
            ResultSetMetaData meta = rs.getMetaData();
            int columnCount = meta.getColumnCount();
            for (int i = 1; i <= columnCount; i++) {
                writer.write(meta.getColumnName(i));
                if (i < columnCount) writer.write("\t");
            }
            writer.newLine();

            // Write rows
            while (rs.next()) {
                for (int i = 1; i <= columnCount; i++) {
                    Object value = rs.getObject(i);
                    writer.write(value == null ? "" : value.toString());
                    if (i < columnCount) writer.write("\t");
                }
                writer.newLine();
                pr.progress();
            }

            Print.printInfo("Deleting parquet " + parquet);
            File parquetFile = new File(parquet);
            parquetFile.delete();
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
}
