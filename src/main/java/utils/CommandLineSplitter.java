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

import java.util.ArrayList;
import java.util.List;

/**
 * Splits one line of command-line options, as a user typed it into a parameter, into the separate
 * arguments a child process is handed.
 *
 * <p>The rules are few on purpose. Quotes and whitespace work much as they do in a shell, and nothing
 * else a shell does is done here, escapes included:
 *
 * <ul>
 *   <li>Whitespace separates arguments.</li>
 *   <li>{@code "..."} or {@code '...'} keeps its text in one argument, whitespace included, and the
 *       quotes themselves are removed. Inside one kind of quote the other kind is an ordinary
 *       character.</li>
 *   <li>There are no escape characters. A backslash is a backslash, because the paths people pass
 *       are Windows paths, and {@code "C:\data\"} has to mean that folder.</li>
 *   <li>Quoted text joins whatever unquoted text it touches, as in a shell:
 *       {@code "C:\My data"\eval.tsv} is the one argument {@code C:\My data\eval.tsv}.</li>
 *   <li>{@code ""} or {@code ''} is an argument of its own that is empty.</li>
 *   <li>A quote that is never closed is an error. Guessing where it was meant to end would hand the
 *       child an argument nobody typed.</li>
 * </ul>
 *
 * <p>Blank text, or none, is no arguments at all, so an empty parameter leaves a command line exactly
 * as it was.
 */
public final class CommandLineSplitter {
    private CommandLineSplitter() {
    }

    /**
     * @param text  the options as typed; {@code null} or blank is no arguments
     * @param param the parameter the text came from, named in the error a bad quote raises
     * @throws IllegalArgumentException when a quote is opened and never closed
     */
    public static List<String> split(String text, String param) {
        List<String> arguments = new ArrayList<>();
        if (text == null) {
            return arguments;
        }
        StringBuilder current = new StringBuilder();
        //an argument can be empty ("" alone) and still be one, so whether one is open is tracked
        //apart from what it holds
        boolean open = false;
        char quote = 0;
        int quotedAt = -1;
        for (int i = 0; i < text.length(); i++) {
            char c = text.charAt(i);
            if (quote != 0) {
                if (c == quote) {
                    quote = 0;
                } else {
                    current.append(c);
                }
            } else if (c == '"' || c == '\'') {
                quote = c;
                quotedAt = i;
                open = true;
            } else if (Character.isWhitespace(c)) {
                if (open) {
                    arguments.add(current.toString());
                    current.setLength(0);
                    open = false;
                }
            } else {
                current.append(c);
                open = true;
            }
        }
        if (quote != 0) {
            throw new IllegalArgumentException(param + " opens a " + quote + " quote at character " +
                    (quotedAt + 1) + " that is never closed. Close it; a quote that belongs to a " +
                    "value goes inside quotes of the other kind, as in \"it's\". " + param + " is: " +
                    text);
        }
        if (open) {
            arguments.add(current.toString());
        }
        return arguments;
    }
}
