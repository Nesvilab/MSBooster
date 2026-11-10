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

package testutils;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.LinkedHashMap;
import java.util.Map;

/** Reads what a Parquet file's footer says, the way a reader of the library would find it. */
public final class ParquetFooters {
    private ParquetFooters() {}

    /** The file's key-value metadata, in file order. */
    public static Map<String, String> read(File parquet) throws SQLException {
        Map<String, String> footer = new LinkedHashMap<>();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             ResultSet rs = stmt.executeQuery("SELECT decode(key) AS k, decode(value) AS v " +
                     "FROM parquet_kv_metadata('" +
                     parquet.getAbsolutePath().replace("\\", "/").replace("'", "''") + "')")) {
            while (rs.next()) {
                footer.put(rs.getString("k"), rs.getString("v"));
            }
        }
        return footer;
    }
}
