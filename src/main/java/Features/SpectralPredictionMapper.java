package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

public interface SpectralPredictionMapper {

    HashMap<String, PredictionEntry> getPreds() throws IOException;
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
                return new mgfFileReader(file, false, executorService, model);
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
