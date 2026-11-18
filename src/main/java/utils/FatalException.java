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

/**
 * A run that cannot go on, for a reason the user can act on, such as a model name MSBooster does not know. The
 * message is written for the user; MainClass prints it as an error, without a stack trace, and exits with status 1.
 * Code on an executor thread throws it like any other exception: the thread that waits on the task gets it back from
 * Future.get, wrapped in an ExecutionException, and {@link #in} finds it there. KoinaException is one for the Koina
 * models.
 */
public class FatalException extends RuntimeException {
    public FatalException(String message) {
        super(message);
    }

    public FatalException(String message, Throwable cause) {
        super(message, cause);
    }

    /** The FatalException that is {@code t} or one of its causes, or null if there is none. */
    public static FatalException in(Throwable t) {
        for (Throwable cause = t; cause != null; cause = cause.getCause()) {
            if (cause instanceof FatalException) {
                return (FatalException) cause;
            }
        }
        return null;
    }
}
