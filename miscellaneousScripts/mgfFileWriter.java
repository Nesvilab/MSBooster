package Features;

public class mgfFileWriter {
    public mgfFileWriter() {

    }

//    public static void writeFineTune(pepXMLReader xml1, pepXMLReader xml2, mzMLReader mzml, String outfile)
//            throws IOException, FileParsingException {
//        fineTuner ft = new fineTuner();
//        ft.addPepXMLReader(xml1);
//        ft.addPepXMLReader(xml2);
//        ft.fineTune();
//
//        FileWriter myWriter = new FileWriter(outfile);
//        for (fineTuner.fineTuneHit hit : ft.hits) {
//            myWriter.write("BEGIN IONS\n");
//            myWriter.write("TITLE=" + hit.pep + "\n");
//            myWriter.write("CHARGE=" + hit.pep.substring(hit.pep.length() - 1) + "\n");
//            //pepmass?
//            myWriter.write("RTINSECONDS=" + hit.rt + "\n");
//
//            //extract mz and intensity
//            double[] MZs = mzml.getMZ(hit.scanNum);
//            double[] Ints = mzml.getIntensity(hit.scanNum);
//            double max_int = (double) Collections.max(Arrays.asList(ArrayUtils.toObject(Ints)));
//            double threshold = max_int / 10.0; //out of memory error
//
//            for (int i = 0; i < MZs.length; i++) {
//                if (Ints[i] > threshold) {
//                    myWriter.write(MZs[i] + "\t" + Ints[i] + "\n");
//                }
//            }
//
//            myWriter.write("END IONS\n");
//        }
//    }

//    public static void main(String[] args) throws IOException, FileParsingException {
//        pepXMLReader x = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/" +
//                "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834_rank1.pepXML");
//        pepXMLReader y = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/" +
//                "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834_rank2.pepXML");
//        mzMLReader mzml = new mzMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/" +
//                "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzML");
//
//        writeFineTune(x, y, mzml, constants.outputFolder + "fineTune.mgf");
//    }
}
