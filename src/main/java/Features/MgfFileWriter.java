package Features;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

public class MgfFileWriter {
    SpectralPredictionMapper spm;
    public MgfFileWriter(SpectralPredictionMapper spm) {
        this.spm = spm;
    }

    public void write(String outfile) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
        for (Map.Entry<String, PredictionEntry> entry : spm.getPreds().entrySet()) {
            String[] peptide = entry.getKey().split("\\|");
            PredictionEntry pe = entry.getValue();
            if (pe.mzs == null) { //was not able to have its mzs predicted
                continue;
            }

            bw.write("BEGIN IONS" + "\n");
            bw.write("TITLE=" + peptide[0] + "\n");
            bw.write("CHARGE=" + peptide[1] + "\n");
            bw.write("RT=" + pe.RT + "\n");
            for (int i = 0; i < pe.mzs.length; i++) {
                if (pe.intensities[i] != 0) {
                    bw.write(pe.mzs[i] + "\t" + pe.intensities[i] + "\n");
                }
            }
            bw.write("END IONS" + "\n");
        }
        bw.close();
    }
}
