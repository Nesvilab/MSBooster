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

package testutils;

import java.security.Permission;

/**
 * Turns a System.exit inside a block into an {@link ExitCalled} exception, so code that exits on an
 * error fails the test that ran it instead of ending the test JVM and every test after it. It needs
 * a JVM that still lets a test install a SecurityManager: Java 17, which prints a deprecation
 * warning for it. On a JVM that refuses, it does nothing and an exit stays an exit.
 */
@SuppressWarnings("removal")
public final class ExitTrap implements AutoCloseable {
    /** Thrown instead of exiting. */
    public static final class ExitCalled extends SecurityException {
        public final int status;

        ExitCalled(int status) {
            super("System.exit(" + status + ") was called");
            this.status = status;
        }
    }

    private final SecurityManager previous;
    private final boolean installed;

    private ExitTrap() {
        previous = System.getSecurityManager();
        boolean ok;
        try {
            System.setSecurityManager(new SecurityManager() {
                @Override
                public void checkPermission(Permission perm) {}

                @Override
                public void checkPermission(Permission perm, Object context) {}

                @Override
                public void checkExit(int status) {
                    throw new ExitCalled(status);
                }
            });
            ok = true;
        } catch (UnsupportedOperationException e) { //Java 18+ without -Djava.security.manager=allow
            ok = false;
        }
        installed = ok;
    }

    public static ExitTrap install() {
        return new ExitTrap();
    }

    @Override
    public void close() {
        if (installed) {
            System.setSecurityManager(previous);
        }
    }
}
