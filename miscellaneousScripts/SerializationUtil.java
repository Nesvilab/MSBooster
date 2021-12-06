package Features;

import umich.ms.fileio.exceptions.FileParsingException;
/**
 * A simple class with generic serialize and deserialize method implementations
 *
 * @author pankaj
 * https://www.journaldev.com/2452/serialization-in-java
 */
public class SerializationUtil {

    // deserialize to Object from given file
    public static Object deserialize(String fileName) throws IOException,
            ClassNotFoundException {
        FileInputStream fis = new FileInputStream(fileName);
        ObjectInputStream ois = new ObjectInputStream(fis);
        Object obj = ois.readObject();
        ois.close();
        return obj;
    }

    // serialize the given object and save it to file
    public static void serialize(Object obj, String fileName)
            throws IOException {
        FileOutputStream fos = new FileOutputStream(fileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(obj);

        fos.close();
    }

    public static void main(String[] args) throws FileParsingException, IOException, ClassNotFoundException {
        long startTime = System.nanoTime();
        mzMLReader m = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/wideWindow/" +
                "23aug2017_hela_serum_timecourse_pool_wide_001.mzML");
        //mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/wideWindowPredsAllRanks.mgf");
        m.createScanNumObjects();
//        System.out.println("loading pepxml");
//        for (int rank = 1; rank < 11; rank ++){
//            pepXMLReader p = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/encyclopedia_params/rank"
//                    + rank + "/23aug2017_hela_serum_timecourse_pool_wide_003.pepXML");
//            m.setPepxmlEntries(p, rank, mgf);
//        }

        try {
            serialize(m, "C:/Users/kevin/OneDriveUmich/proteomics/wide1SerializeTest");
        } catch (IOException e){
            e.printStackTrace();
            return;
        }

//        mzMLReader m = (mzMLReader) deserialize("C:/Users/kevin/OneDriveUmich/proteomics/wide2Serialize");
//        m.setScans();
//        System.out.println("loading pepxml");
//        for (int rank = 1; rank < 11; rank ++){
//            pepXMLReader p = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/encyclopedia_params/rank"
//                    + rank + "/23aug2017_hela_serum_timecourse_pool_wide_001.pepXML");
//            m.setPepxmlEntries(p, rank, mgf);
//        }
//
//        long endTime = System.nanoTime();
//        long duration = (endTime - startTime) / 1000000;
//        System.out.println(duration);
    }
}
