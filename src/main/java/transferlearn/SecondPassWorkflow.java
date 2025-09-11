package transferlearn;

import utils.Print;

public class SecondPassWorkflow {
    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.SecondPassWorkflow " +
                "--paramsList <msbooster parameters> " +
                "--urlTransferLearn <server url with transfer learning port> --urlPredict <server url with predict port> " +
                "--library <path/to/librarytsv> " +
                "optional: --ms2 <predict ms2> --rt <predict rt> --im <predict ccs> --basename <output base name>");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.SecondPassWorkflow " +
                "--paramsList msbooster_params.txt " +
                "--urlTransferLearn http://localhost:8000 --urlPredict http://localhost:8001 " +
                "--library library.tsv --ms2 true --rt true --im false");
    }

    public static void main(String[] args) throws Exception {
        //parse arguments
        if (args.length % 2 != 0) {
            errorMessage();
        }

        String params = "";
        String urlTransferLearn = "";
        String urlPredict = "";
        String library = "";
        String ms2 = "true";
        String rt = "true";
        String im = "true";
        String basename = "";

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--paramsList":
                    params =  args[i + 1];
                    break;
                case "--urlTransferLearn":
                    urlTransferLearn =  args[i + 1];
                    break;
                case "--urlPredict":
                    urlPredict =  args[i + 1];
                    break;
                case "--library":
                    library =  args[i + 1];
                    break;
                case "--ms2":
                    ms2 = args[i + 1];
                    break;
                case "--rt":
                    rt = args[i + 1];
                    break;
                case "--im":
                    im = args[i + 1];
                    break;
                case "--basename":
                    basename =  args[i + 1];
                    break;
            }
        }

        if (params.isEmpty() || urlTransferLearn.isEmpty() || urlPredict.isEmpty() || library.isEmpty()) {
            errorMessage();
        }

        String modelZip = Trainer.train(new String[] {
                "--url", urlTransferLearn,
                "--library", library,
                "--basename", basename});
        Predictor.main(new String[] {
                "--paramsList", params,
                "--url", urlPredict,
                "--model", modelZip,
                "--ms2", ms2,
                "--rt", rt,
                "--im", im
        });
    }
}
