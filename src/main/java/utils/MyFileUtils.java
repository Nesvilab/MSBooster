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

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class MyFileUtils {
    public static void deleteWholeDirectory(String directory) {
        try {
            Path rootPath = Paths.get(directory);
            if (Files.exists(rootPath)) {
                try (Stream<Path> walk = Files.walk(rootPath)) {
                    walk.sorted(Comparator.reverseOrder())
                            .map(Path::toFile)
                            .forEach(File::delete);
                }
            }
        } catch (Exception ignored) {}
    }

    public static void createWholeDirectory(String directory) {
        if (! new File(directory).exists()) {
            new File(directory).mkdirs();
        }
    }

    //solution from chatgpt
    public static String findDeepestCommonDirectory(String[] filePaths) {
        Path root = Paths.get(filePaths[0]).toAbsolutePath().getRoot();

        List<List<String>> splitPaths = new ArrayList<>();

        for (String pathStr : filePaths) {
            Path parentPath = Paths.get(pathStr).toAbsolutePath().normalize().getParent();
            List<String> segments = new ArrayList<>();
            for (Path part : parentPath) {
                segments.add(part.toString());
            }
            splitPaths.add(segments);
        }

        List<String> commonSegments = new ArrayList<>();

        for (int i = 0; ; i++) {
            String segment = null;
            for (List<String> segments : splitPaths) {
                if (i >= segments.size()) {
                    return buildPathString(commonSegments, root);
                }
                if (segment == null) {
                    segment = segments.get(i);
                } else if (!segment.equals(segments.get(i))) {
                    return buildPathString(commonSegments, root);
                }
            }
            commonSegments.add(segment);
        }
    }

    private static String buildPathString(List<String> segments, Path root) {
        Path fullPath = root;
        for (String segment : segments) {
            fullPath = fullPath.resolve(segment);
        }
        return fullPath.toString();
    }
}
