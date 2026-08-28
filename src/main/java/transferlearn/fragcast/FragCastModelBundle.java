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

package transferlearn.fragcast;

import allconstants.Constants;
import allconstants.FragCastWeights;
import utils.Print;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

/**
 * A fine-tuned FragCast model as one file.
 *
 * <p>FragCast keeps a separate ONNX file per property - one for retention time, one for ion
 * mobility, one for MS2 - and takes each as its own command-line flag. That is the right shape for
 * the executable and the wrong shape for a person: three files have to be kept together, kept
 * straight, and named individually every time the model is used. This class is the seam between the
 * two. A fine-tune writes one zip; a later run opens it and gets the three paths back.
 *
 * <p>The zip is flat and self-describing:
 *
 * <pre>
 *   manifest.properties
 *   FragCast-RT.onnx     (optional)
 *   FragCast-IM.onnx     (optional)
 *   FragCast-Spec.onnx   (optional)
 * </pre>
 *
 * <p>The manifest is not decoration. The RT and IM models are the same architecture and the same
 * shape, so a file swapped for the other loads without complaint and quietly predicts the wrong
 * property; nothing in the ONNX itself says which is which. The entry name is what says it, and the
 * checksum beside it is what says the file is the one that name means - a zip does not enforce its
 * own CRC on read, so a corrupted entry inflates to different bytes with no error at all.
 *
 * <p>Nothing here records a path. A bundle names its own entries and no directory, so it stays valid
 * wherever it is copied to.
 */
public final class FragCastModelBundle {

    /** Bumped only for a change an older MSBooster could not read; it refuses anything higher. */
    public static final int BUNDLE_VERSION = 1;

    public static final String MANIFEST_NAME = "manifest.properties";

    /** Bundle entry names, one per task. The name inside the bundle never varies. */
    private static final Map<String, String> ENTRY_NAMES = new LinkedHashMap<>();

    static {
        ENTRY_NAMES.put("rt", "FragCast-RT.onnx");
        ENTRY_NAMES.put("im", "FragCast-IM.onnx");
        ENTRY_NAMES.put("spec", "FragCast-Spec.onnx");
    }

    /** A bundle holds at most one manifest and one model per task. */
    private static final int MAX_ENTRIES = FragCastWeights.TASKS.length + 1;

    /** Comfortably above the largest real model (2.9 MB) and far below anything alarming. */
    private static final long MAX_ENTRY_BYTES = 64L * 1024 * 1024;

    private FragCastModelBundle() {
    }

    // ---------------------------------------------------------------------------------------
    // writing
    // ---------------------------------------------------------------------------------------

    /**
     * Write every model {@code weights} names into one zip.
     *
     * <p>Every named model goes in, not only the ones this run fine-tuned. A bundle that carried
     * just the fine-tuned properties would predict differently from the run that produced it as soon
     * as one of the others had been supplied by hand, and the point of the file is that it
     * reproduces that run.
     *
     * @param weights the models to bundle; empty properties are skipped
     * @param dest    the zip to write, overwritten if it exists
     * @return {@code dest}, or {@code null} when {@code weights} names no model at all
     */
    public static File write(FragCastWeights weights, File dest) throws IOException {
        if (weights == null || weights.isBase()) {
            return null;
        }
        final Map<String, File> models = new LinkedHashMap<>();
        for (String task : FragCastWeights.TASKS) {
            final String path = weights.forTask(task);
            if (path.isEmpty()) {
                continue;
            }
            final File file = new File(path);
            if (!file.isFile()) {
                throw new IOException("cannot bundle the FragCast " + task + " model: no such file " + path);
            }
            models.put(task, file);
        }
        final Map<String, String> manifest = new LinkedHashMap<>();
        manifest.put("bundleVersion", String.valueOf(BUNDLE_VERSION));
        manifest.put("producer", Constants.versionNumber);
        manifest.put("createdUtc", Instant.now().toString());
        manifest.put("tasks", String.join(",", models.keySet()));

        final File parent = dest.getAbsoluteFile().getParentFile();
        if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
            throw new IOException("could not create " + parent);
        }

        try (OutputStream out = Files.newOutputStream(dest.toPath());
             ZipOutputStream zip = new ZipOutputStream(out)) {
            for (Map.Entry<String, File> e : models.entrySet()) {
                final String task = e.getKey();
                final File file = e.getValue();
                final byte[] bytes = Files.readAllBytes(file.toPath());
                final String entryName = ENTRY_NAMES.get(task);

                manifest.put(task + ".sha256", sha256(bytes));

                zip.putNextEntry(new ZipEntry(entryName));
                zip.write(bytes);
                zip.closeEntry();
            }
            // last, so a reader that streams the zip has seen every file the manifest describes
            zip.putNextEntry(new ZipEntry(MANIFEST_NAME));
            final Properties props = new Properties();
            props.putAll(manifest);
            final Writer manifestWriter = new OutputStreamWriter(zip, StandardCharsets.UTF_8);
            props.store(manifestWriter, "FragCast model bundle. See transferlearn.fragcast.FragCastModelBundle.");
            manifestWriter.flush();
            zip.closeEntry();
        }
        return dest;
    }

    // ---------------------------------------------------------------------------------------
    // reading
    // ---------------------------------------------------------------------------------------

    /**
     * Extract {@code zip} into {@code destDir} and return the models it holds.
     *
     * <p>Everything is checked before anything is believed. A bundle that fails any check is refused
     * whole rather than in part: a truncated or edited archive is not a model with one bad file in
     * it, it is an archive nobody should predict from.
     *
     * @throws IOException with a message naming what was wrong, ready to print at the user
     */
    public static FragCastWeights open(File zip, File destDir) throws IOException {
        if (zip == null || !zip.isFile()) {
            throw new IOException("no such FragCast model file: " + zip);
        }
        if (!destDir.isDirectory() && !destDir.mkdirs()) {
            throw new IOException("could not create " + destDir + " to unpack " + zip.getName() + " into");
        }

        try (ZipFile archive = openZip(zip)) {
            final Map<String, ZipEntry> entries = readEntries(zip, archive);
            final Properties manifest = readManifest(zip, archive, entries);
            final List<String> tasks = declaredTasks(zip, manifest);

            final Set<String> referenced = new LinkedHashSet<>();
            referenced.add(MANIFEST_NAME);
            for (String task : tasks) {
                referenced.add(ENTRY_NAMES.get(task));
            }
            for (String name : entries.keySet()) {
                if (!referenced.contains(name)) {
                    // not pedantry: an unreferenced file is either a bundle built by something that
                    // disagrees with this reader, or a passenger, and neither should be unpacked
                    throw new IOException(describe(zip) + " holds " + name + ", which its manifest does not mention");
                }
            }

            FragCastWeights weights = FragCastWeights.base();
            for (String task : tasks) {
                final String entryName = ENTRY_NAMES.get(task);
                final ZipEntry entry = entries.get(entryName);
                if (entry == null) {
                    throw new IOException(describe(zip) + " names " + entryName + " for its " + task +
                            " model but does not contain it");
                }
                final byte[] bytes = read(zip, archive, entry);
                verify(zip, manifest, task, bytes);

                final File out = new File(destDir, entryName);
                Files.write(out.toPath(), bytes);
                weights = weights.withTask(task, out.getAbsolutePath());
            }
            return weights;
        }
    }

    private static ZipFile openZip(File zip) throws IOException {
        try {
            return new ZipFile(zip);
        } catch (IOException e) {
            throw new IOException(zip + " is not a readable zip file (" + e.getMessage() + ")");
        }
    }

    /** Entry names, checked to be flat and sane before any of them is opened. */
    private static Map<String, ZipEntry> readEntries(File zip, ZipFile archive) throws IOException {
        final Map<String, ZipEntry> entries = new LinkedHashMap<>();
        final Enumeration<? extends ZipEntry> it = archive.entries();
        while (it.hasMoreElements()) {
            final ZipEntry entry = it.nextElement();
            final String name = entry.getName();
            if (entry.isDirectory() || name.indexOf('/') >= 0 || name.indexOf('\\') >= 0
                    || name.contains("..") || name.trim().isEmpty()) {
                // a bundle is flat by construction, so anything that could escape destDir is not a
                // path to be sanitised - it is a bundle this reader did not write
                throw new IOException(describe(zip) + " holds the file " + name +
                        ", but a FragCast model bundle contains only plain file names");
            }
            if (entries.put(name, entry) != null) {
                throw new IOException(describe(zip) + " holds two files named " + name);
            }
            if (entries.size() > MAX_ENTRIES) {
                throw new IOException(describe(zip) + " holds more than " + MAX_ENTRIES + " files");
            }
            final long size = entry.getSize();
            if (size > MAX_ENTRY_BYTES) {
                throw new IOException(describe(zip) + " holds " + name + " at " + size +
                        " bytes, past the " + MAX_ENTRY_BYTES + " byte limit for a model file");
            }
        }
        return entries;
    }

    private static Properties readManifest(File zip, ZipFile archive, Map<String, ZipEntry> entries)
            throws IOException {
        final ZipEntry entry = entries.get(MANIFEST_NAME);
        if (entry == null) {
            throw new IOException(describe(zip) + " has no " + MANIFEST_NAME +
                    ", so it is not a FragCast model bundle");
        }
        final Properties manifest = new Properties();
        try (InputStream in = archive.getInputStream(entry)) {
            manifest.load(in);
        }
        final String declared = manifest.getProperty("bundleVersion", "").trim();
        final int version;
        try {
            version = Integer.parseInt(declared);
        } catch (NumberFormatException e) {
            throw new IOException(describe(zip) + " declares the unreadable bundle version '" + declared + "'");
        }
        if (version > BUNDLE_VERSION) {
            throw new IOException(describe(zip) + " is a version " + version + " bundle, and this is " +
                    Constants.versionNumber + ", which reads up to version " + BUNDLE_VERSION +
                    ". Use the MSBooster version that wrote it.");
        }
        return manifest;
    }

    private static List<String> declaredTasks(File zip, Properties manifest) throws IOException {
        final List<String> tasks = new ArrayList<>();
        for (String raw : requireValue(zip, manifest, "tasks").split(",")) {
            final String task = raw.trim().toLowerCase();
            if (task.isEmpty()) {
                continue;
            }
            if (!ENTRY_NAMES.containsKey(task)) {
                throw new IOException(describe(zip) + " names the unknown task '" + task + "'");
            }
            if (tasks.contains(task)) {
                throw new IOException(describe(zip) + " names the task '" + task + "' twice");
            }
            tasks.add(task);
        }
        if (tasks.isEmpty()) {
            throw new IOException(describe(zip) + " holds no model");
        }
        return tasks;
    }

    /** Read one entry, refusing to grow past the size its own header declared. */
    private static byte[] read(File zip, ZipFile archive, ZipEntry entry) throws IOException {
        try (InputStream in = archive.getInputStream(entry)) {
            final java.io.ByteArrayOutputStream buffer = new java.io.ByteArrayOutputStream();
            final byte[] chunk = new byte[8192];
            long total = 0;
            int n;
            while ((n = in.read(chunk)) > 0) {
                total += n;
                if (total > MAX_ENTRY_BYTES) {
                    throw new IOException(describe(zip) + " expands " + entry.getName() +
                            " past the " + MAX_ENTRY_BYTES + " byte limit for a model file");
                }
                buffer.write(chunk, 0, n);
            }
            return buffer.toByteArray();
        }
    }

    private static void verify(File zip, Properties manifest, String task, byte[] bytes) throws IOException {
        final String name = ENTRY_NAMES.get(task);

        //no length check: the digest is taken over the bytes, so a truncated entry fails it too, and
        //parsing a length out of the manifest let an unchecked NumberFormatException escape open()
        final String declaredHash = requireValue(zip, manifest, task + ".sha256").trim().toLowerCase();
        final String actualHash = sha256(bytes);
        if (!declaredHash.equals(actualHash)) {
            throw new IOException(describe(zip) + " declares " + name + " with checksum " + declaredHash +
                    " but holds " + actualHash + "; the bundle is damaged or was edited");
        }
    }

    // ---------------------------------------------------------------------------------------
    // wiring a bundle into a run
    // ---------------------------------------------------------------------------------------

    /**
     * Unpack {@link Constants#FragCastModelZip}, if one was named, and point the run's FragCast
     * weights at what it holds. Does nothing when no bundle was named.
     *
     * <p>A property the parameters already name explicitly wins over the bundle, and the run says so
     * rather than overriding it silently: someone who set {@code FragCastSpecOnnx} by hand alongside
     * a bundle meant that file, and a bundle that quietly replaced it would be very hard to notice.
     *
     * @throws IOException if the bundle cannot be read. Naming a bundle is explicit, so an unusable
     *                     one is an error rather than something to carry on past.
     */
    public static void applyToConstants() throws IOException {
        final String zipPath = Constants.FragCastModelZip == null ? "" : Constants.FragCastModelZip.trim();
        if (zipPath.isEmpty() || zipPath.equals("null")) {
            return;
        }
        final File zip = new File(zipPath);
        final File destDir = unpackDirectory(zip);
        Print.printInfo("Opening the FragCast model " + zip.getName() + " in " + destDir);

        final FragCastWeights bundled = open(zip, destDir);
        final FragCastWeights named = FragCastWeights.fromConstants();

        for (String task : FragCastWeights.TASKS) {
            final String fromBundle = bundled.forTask(task);
            if (fromBundle.isEmpty()) {
                continue;
            }
            final String alreadyNamed = named.forTask(task);
            if (!alreadyNamed.isEmpty()) {
                Print.printInfo("  keeping the " + task + " model named by " + paramName(task) + " = " +
                        alreadyNamed + "; the one in " + zip.getName() + " is not used");
                continue;
            }
            setConstant(task, fromBundle);
            Print.printInfo("  " + task + " model: " + fromBundle);
        }
    }

    /**
     * Where a bundle is unpacked: one fixed directory inside the run's own output, never a fresh
     * temporary one.
     *
     * <p>The paths matter beyond finding the files. {@link FragCastWeights#fileTag()} hashes them
     * into the name of the prediction file the run writes, so a directory that changed per
     * invocation would change that name too and every rerun would re-predict a library it already
     * had.
     */
    static File unpackDirectory(File zip) {
        final String output = Constants.outputDirectory;
        if (output != null && !output.trim().isEmpty()) {
            return new File(output, "fragcast_model");
        }
        // no output directory yet - the library-only workflows never establish one - so the bundle's
        // own directory is the one location the caller has definitely named
        final File parent = zip.getAbsoluteFile().getParentFile();
        return new File(parent, "fragcast_model");
    }

    // ---------------------------------------------------------------------------------------
    // helpers
    // ---------------------------------------------------------------------------------------

    private static String sha256(byte[] bytes) {
        final MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is required of every Java runtime", e);
        }
        final StringBuilder hex = new StringBuilder(64);
        for (byte b : digest.digest(bytes)) {
            hex.append(Character.forDigit((b >> 4) & 0xf, 16));
            hex.append(Character.forDigit(b & 0xf, 16));
        }
        return hex.toString();
    }

    private static String requireValue(File zip, Properties manifest, String key) throws IOException {
        final String value = manifest.getProperty(key, "").trim();
        if (value.isEmpty()) {
            throw new IOException(describe(zip) + " is missing '" + key + "' from its " + MANIFEST_NAME);
        }
        return value;
    }

    private static String describe(File zip) {
        return "The FragCast model " + zip.getName();
    }

    private static void setConstant(String task, String path) {
        switch (task) {
            case "rt":
                Constants.FragCastRtOnnx = path;
                break;
            case "im":
                Constants.FragCastImOnnx = path;
                break;
            case "spec":
                Constants.FragCastSpecOnnx = path;
                break;
            default:
                break;
        }
    }

    private static String paramName(String task) {
        switch (task) {
            case "rt":
                return "FragCastRtOnnx";
            case "im":
                return "FragCastImOnnx";
            default:
                return "FragCastSpecOnnx";
        }
    }
}
