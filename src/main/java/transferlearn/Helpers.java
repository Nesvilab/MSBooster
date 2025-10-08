package transferlearn;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import utils.Print;

import javax.net.ssl.HttpsURLConnection;
import java.io.*;
import java.lang.reflect.Type;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.sql.*;
import java.util.HashMap;

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

    static void convertCsvToParquet(String csvFilePath, String parquetFilePath) {
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

        } catch (SQLException e) {
            Print.printError("Error converting " + csvFilePath + " to parquet: " + e.getMessage());
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    static void convertParquetToCsv(String parquet, String tsv) throws SQLException, IOException {
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             BufferedWriter writer = new BufferedWriter(new FileWriter(tsv))) {

            // Register and query Parquet
            String query = "SELECT * FROM read_parquet('" + parquet + "')";
            ResultSet rs = stmt.executeQuery(query);

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
            }
        }
    }
}
