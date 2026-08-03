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
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class CaseInsensitiveHashSet extends HashSet<String> {
    public CaseInsensitiveHashSet(String[] elements) {
        super(Arrays.asList(elements));
    }
    public CaseInsensitiveHashSet(List<String> elements) {
        super(elements);
    }
    public CaseInsensitiveHashSet() {};

    @Override
    public boolean contains(Object o) {
        if (!(o instanceof String))
            return false;
        String element = (String) o;
        for (String s : this) {
            if (s.equalsIgnoreCase(element))
                return true;
        }
        return false;
    }
}

