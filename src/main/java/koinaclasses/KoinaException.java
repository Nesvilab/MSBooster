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

package koinaclasses;

import utils.FatalException;

/**
 * A Koina prediction that cannot go on. The message is written for the user and names the model; MainClass prints
 * it and exits with status 1, as for any {@link FatalException}. Code on an executor thread throws it like any other
 * exception: the thread that waits on the task gets it back from Future.get, wrapped in an ExecutionException, and
 * {@link #in} finds it there.
 */
public class KoinaException extends FatalException {
    public KoinaException(String message) {
        super(message);
    }

    public KoinaException(String message, Throwable cause) {
        super(message, cause);
    }

    /** The KoinaException that is {@code t} or one of its causes, or null if there is none. */
    public static KoinaException in(Throwable t) {
        for (Throwable cause = t; cause != null; cause = cause.getCause()) {
            if (cause instanceof KoinaException) {
                return (KoinaException) cause;
            }
        }
        return null;
    }
}
