package utils;

import java.util.Arrays;
import java.util.HashSet;

public class CaseInsensitiveHashSet extends HashSet<String> {
    public CaseInsensitiveHashSet(String[] elements) {
        super(Arrays.asList(elements));
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

