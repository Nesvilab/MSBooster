/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
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

