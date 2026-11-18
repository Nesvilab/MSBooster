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

package mainsteps;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import allconstants.LowercaseModelMapper;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import koinaclasses.KoinaException;
import koinaclasses.KoinaModelSpec.Kind;
import org.junit.jupiter.api.Test;
import utils.FatalException;
import writers.JSONWriter;

/**
 * MainClass ends a run that threw with exit status 1. A FatalException, such as a KoinaException, also one thrown
 * on an executor thread and handed back wrapped in an ExecutionException, is a message for the user, printed as an
 * error without a stack trace.
 */
class MainClassFatalTest {
    //a model name MSBooster does not know, such as the PredFull it no longer has, ends the run with status 1 (it
    //ended it with status 0 before)
    @Test
    void anUnknownModelIsFatal() {
        Map<String, String> names = new LowercaseModelMapper().getLowercaseToModel();
        FatalException e = assertThrows(FatalException.class,
                () -> MainClass.modelName("spectra", "spectraModel", "PredFull", names, Kind.MS2));
        assertTrue(e.getMessage().startsWith("No spectra model called PredFull (spectraModel). The Koina spectra "
                + "models MSBooster can call are "), e.getMessage());
        assertTrue(e.getMessage().contains("Prosit_2020_intensity_HCD"), e.getMessage());
        assertFalse(e.getMessage().contains("Deeplc_hela_hf"), "only the models of the parameter's kind");

        String[] printed = printFatal(e);
        assertTrue(printed[0].contains("[ERROR] - No spectra model called PredFull"), printed[0]);
        assertTrue(printed[0].contains("[ERROR] - Exiting"), printed[0]);
        assertEquals("", printed[1], "no stack trace for a message to the user");

        assertThrows(FatalException.class,
                () -> MainClass.modelName("RT", "rtModel", "Deeplc", names, Kind.RT));
        assertThrows(FatalException.class,
                () -> MainClass.modelName("IM", "imModel", "PredFull", names, Kind.IM));
    }

    @Test
    void aKnownModelIsSpelledAsMSBoosterSpellsIt() {
        Map<String, String> names = new LowercaseModelMapper().getLowercaseToModel();
        assertEquals("UniSpec", MainClass.modelName("spectra", "spectraModel", "unispec", names, Kind.MS2));
        assertEquals("Deeplc_hela_hf", MainClass.modelName("RT", "rtModel", "DEEPLC_HELA_HF", names, Kind.RT));
        assertEquals("DIA-NN", MainClass.modelName("IM", "imModel", "dia-nn", names, Kind.IM));
        assertEquals("", MainClass.modelName("auxiliary spectra", "auxSpectraModel", "", names, Kind.MS2),
                "no model");
    }

    //the request writer is only reached with a Koina model; any other name is a KoinaException, not an exit
    @Test
    void theRequestWriterRefusesAModelKoinaDoesNotHave() {
        String savedAux = Constants.auxSpectraModel;
        Constants.auxSpectraModel = "";
        try {
            KoinaException e = assertThrows(KoinaException.class,
                    () -> new JSONWriter("PredFull", Set.of("PEPTIDEK,2,30,QE,HCD"), false, 100));
            assertEquals("PredFull is not a Koina model MSBooster can call", e.getMessage());
        } finally {
            Constants.auxSpectraModel = savedAux;
        }
    }

    @Test
    void aKoinaExceptionIsPrintedAsAnError() {
        KoinaException koina = new KoinaException("MSBooster could not read the predictions of Koina model X");
        String[] printed = printFatal(new ExecutionException(new RuntimeException("task", koina)));

        assertTrue(printed[0].contains("[ERROR] - MSBooster could not read the predictions of Koina model X"),
                printed[0]);
        assertTrue(printed[0].contains("[ERROR] - Exiting"), printed[0]);
        assertEquals("", printed[1], "no stack trace for a message to the user");
    }

    @Test
    void anyOtherExceptionKeepsItsStackTrace() {
        String[] printed = printFatal(new IllegalStateException("a bug"));
        assertFalse(printed[0].contains("[ERROR]"), printed[0]);
        assertTrue(printed[1].contains("java.lang.IllegalStateException: a bug"), printed[1]);
        assertTrue(printed[1].contains("\tat "), printed[1]);
    }

    @Test
    void theKoinaExceptionIsFoundAmongTheCauses() {
        KoinaException koina = new KoinaException("inner");
        assertEquals(koina, KoinaException.in(new ExecutionException(new RuntimeException(koina))));
        assertEquals(koina, KoinaException.in(koina));
        assertEquals(null, KoinaException.in(new RuntimeException(new IllegalStateException())));
    }

    //what MainClass.printFatal writes to System.out and to System.err
    private static String[] printFatal(Throwable e) {
        PrintStream out = System.out;
        PrintStream err = System.err;
        ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream errBytes = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(errBytes, true, StandardCharsets.UTF_8));
            MainClass.printFatal(e);
        } finally {
            System.setOut(out);
            System.setErr(err);
        }
        return new String[] {outBytes.toString(StandardCharsets.UTF_8), errBytes.toString(StandardCharsets.UTF_8)};
    }
}
