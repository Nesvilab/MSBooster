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

package utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

// What a user types into a FragCast options box reaches FragCast as these arguments and nothing
// else: a token split in the wrong place is a flag FragCast refuses, or worse, a value it accepts
// for the wrong flag.
public class CommandLineSplitterTest {
    private static List<String> split(String text) {
        return CommandLineSplitter.split(text, "FragCastFineTuneCmdOpts");
    }

    @Test
    public void whitespaceOfAnyKindAndAmountSeparatesArguments() {
        assertEquals(List.of("--epochs", "20", "--lr", "1e-4"),
                split("  --epochs 20\t--lr   1e-4 "));
    }

    // An empty parameter must leave the command exactly as it was, not add an empty argument.
    @Test
    public void blankTextIsNoArguments() {
        assertEquals(List.of(), split(""));
        assertEquals(List.of(), split("   \t "));
        assertEquals(List.of(), split(null));
    }

    @Test
    public void eitherQuoteKeepsSpacesInOneArgumentAndIsRemoved() {
        assertEquals(List.of("--eval", "C:\\My data\\eval.tsv"),
                split("--eval \"C:\\My data\\eval.tsv\""));
        assertEquals(List.of("--eval", "C:\\My data\\eval.tsv"),
                split("--eval 'C:\\My data\\eval.tsv'"));
    }

    // Windows paths are what people pass, so a backslash escapes nothing - not even a quote, which
    // is what lets a quoted folder end in one.
    @Test
    public void aBackslashIsJustABackslash() {
        assertEquals(List.of("--eval", "C:\\data\\", "--x", "a\\b"),
                split("--eval \"C:\\data\\\" --x a\\b"));
    }

    @Test
    public void quotedTextJoinsTheUnquotedTextItTouches() {
        assertEquals(List.of("C:\\My data\\eval.tsv"), split("\"C:\\My data\"\\eval.tsv"));
        assertEquals(List.of("abc d"), split("a\"b\"'c d'"));
    }

    @Test
    public void theOtherKindOfQuoteIsAnOrdinaryCharacterInsideOne() {
        assertEquals(List.of("it's", "say \"hi\""), split("\"it's\" 'say \"hi\"'"));
    }

    // "" is how a user hands FragCast an empty value; dropping it would shift every argument after it.
    @Test
    public void anEmptyPairOfQuotesIsAnEmptyArgument() {
        assertEquals(List.of("--decoy-prefix", "", "--epochs", "5"),
                split("--decoy-prefix \"\" --epochs 5"));
        assertEquals(List.of(""), split("''"));
    }

    // Guessing where an unclosed quote was meant to end would hand FragCast an argument nobody typed,
    // so it is refused, naming the parameter the user has to go and fix.
    @Test
    public void anUnterminatedQuoteIsAnErrorNamingTheParameter() {
        for (String text : new String[] {"--rt-onnx \"C:\\My models\\RT.onnx", "--epochs 5 'x", "\""}) {
            IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                    () -> CommandLineSplitter.split(text, "FragCastLibraryCmdOpts"));
            assertTrue(e.getMessage().startsWith("FragCastLibraryCmdOpts "), e.getMessage());
            assertTrue(e.getMessage().contains("never closed"), e.getMessage());
            assertTrue(e.getMessage().endsWith(text), "the text was not quoted back: " + e.getMessage());
        }
    }
}
