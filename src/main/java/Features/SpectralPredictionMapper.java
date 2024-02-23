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

package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

public interface SpectralPredictionMapper {

    ConcurrentHashMap<String, PredictionEntry> getPreds() throws IOException;
    void clear();

    float getMaxPredRT();

    //TODO: multithread all of these
    static SpectralPredictionMapper createSpectralPredictionMapper(String file, String model,
                                                                   ExecutorService executorService)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {
        //detecting file extension
        String[] extensionSplit = file.split("\\.");
        String extension = extensionSplit[extensionSplit.length - 1];

        switch (extension) {
            case "bin":
                return new DiannSpeclibReader(file);
            case "mgf":
                //for pdeep and alphapeptdeep
                return new MgfFileReader(file, false, executorService, model);
            case "msp":
                return new MspReader(file);
            case "dlib":
                return new DlibReader(file);
            default:
                System.out.println(extension + " is not a valid prediction file format");
                return null;
        }
    }

    static SpectralPredictionMapper createSpectralPredictionMapper(String file, File[] pinFiles,
                                                                   ExecutorService executorService)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {
        return new PredFullSpeclibReader(file, false, pinFiles, executorService);
    }
}
