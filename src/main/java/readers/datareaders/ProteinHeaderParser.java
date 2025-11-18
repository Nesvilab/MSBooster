package readers.datareaders;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ProteinHeaderParser {
    
    public enum DatabaseType {
        UNIPROT, NCBI, CPTAC_ENSEMBL, ENSEMBL, UNIREF, TAIR, NEXTPROT, GENERIC
    }
    
    public static class ProteinRecord {
        public String originalHeader;
        public String partHeader;
        public String id;
        String entryName;
        public String geneName;
        String proteinName;
        String organism;
        public DatabaseType dbType;
        public boolean isDecoy;
    }
    
    private static final Pattern GENE_NAME_ENSEMBL = Pattern.compile("(ENSG\\d{1,11}\\.?\\d?\\d?)");
    private static final Pattern GENE_NAME_CPTAC_ENSEMBL = Pattern.compile("(ENSG\\d{1,11}\\.?\\d?\\d?)");
    private static final Pattern GENE_NAME_NCBI = Pattern.compile("GN\\=(.+)\\s[\\[|OX\\=|GN\\=|PE\\=|$]");
    private static final Pattern GENE_NAME_UNIPROT = Pattern.compile("GN\\=([\\S]+)");
    private static final Pattern GENE_NAME_TAIR = Pattern.compile("\\|\\sSymbols:\\s*(.+?)\\s\\|");
    
    private static final Pattern ID_UNIPROT = Pattern.compile("[sp|tr]\\|(.+?)\\|");
    private static final Pattern ID_NCBI = Pattern.compile("(\\w{2}_\\d{1,10}\\.?(\\d{1,2})?)");
    private static final Pattern ID_ENSEMBL = Pattern.compile("(ENSP\\w+(?:\\.\\d+)?)");
    private static final Pattern ID_NEXTPROT = Pattern.compile("nxp\\|(.+?)\\|");
    
    public static DatabaseType classify(String header, String decoyTag) {
        String seq = header;
        if (seq.startsWith(decoyTag)) {
            seq = seq.substring(decoyTag.length());
        }
        if (seq.startsWith("contam_")) {
            seq = seq.substring(7);
        }
        if (seq.startsWith(decoyTag)) {
            seq = seq.substring(decoyTag.length());
        }
        
        if (seq.startsWith("sp|") || seq.startsWith("tr|") || seq.startsWith("db|")) {
            return DatabaseType.UNIPROT;
        } else if (seq.startsWith("AP_") || seq.startsWith("NP_") || 
                   seq.startsWith("YP_") || seq.startsWith("XP_") || 
                   seq.startsWith("ZP") || seq.startsWith("WP_")) {
            return DatabaseType.NCBI;
        } else if (seq.contains("ENSP") && seq.contains("|ENST") && seq.contains("|ENSG")) {
            return DatabaseType.CPTAC_ENSEMBL;
        } else if (seq.startsWith("ENSP")) {
            return DatabaseType.ENSEMBL;
        } else if (seq.startsWith("UniRef")) {
            return DatabaseType.UNIREF;
        } else if (seq.startsWith("AT")) {
            if (seq.matches("^AT.+\\s\\|\\sSymbols:[^\\|]+\\s\\|\\s[^\\|]+\\s\\|.*")) {
                return DatabaseType.TAIR;
            }
        } else if (seq.startsWith("nxp")) {
            return DatabaseType.NEXTPROT;
        }
        
        return DatabaseType.GENERIC;
    }
    
    public static String getGeneName(String header, DatabaseType dbType) {
        Pattern pattern = null;
        
        switch (dbType) {
            case ENSEMBL:
                pattern = GENE_NAME_ENSEMBL;
                break;
            case CPTAC_ENSEMBL:
                pattern = GENE_NAME_CPTAC_ENSEMBL;
                break;
            case NCBI:
                pattern = GENE_NAME_NCBI;
                break;
            case UNIPROT:
                pattern = GENE_NAME_UNIPROT;
                break;
            case TAIR:
                pattern = GENE_NAME_TAIR;
                break;
            case NEXTPROT:
                String[] parts = header.split("\\|");
                if (parts.length > 2) {
                    return parts[2].trim();
                }
                return "";
            case UNIREF:
            case GENERIC:
            default:
                return "";
        }
        
        Matcher matcher = pattern.matcher(header);
        if (matcher.find() && matcher.groupCount() >= 1) {
            return matcher.group(1).trim();
        }
        
        return "";
    }
    
    public static String getID(String header, DatabaseType dbType) {
        Pattern pattern = null;
        
        switch (dbType) {
            case UNIPROT:
                pattern = ID_UNIPROT;
                break;
            case NCBI:
                pattern = ID_NCBI;
                break;
            case ENSEMBL:
            case CPTAC_ENSEMBL:
                pattern = ID_ENSEMBL;
                break;
            case NEXTPROT:
                pattern = ID_NEXTPROT;
                break;
            case GENERIC:
                return header;
            default:
                return "";
        }
        
        Matcher matcher = pattern.matcher(header);
        if (matcher.find() && matcher.groupCount() >= 1) {
            return matcher.group(1).trim();
        }
        
        return header;
    }
    
    public static ProteinRecord processHeader(String header, String decoyTag) {
        ProteinRecord record = new ProteinRecord();
        
        record.originalHeader = header;
        int idx = header.indexOf(" ");
        record.partHeader = (idx == -1) ? header : header.substring(0, idx);
        
        record.isDecoy = header.startsWith(decoyTag);
        record.dbType = classify(header, decoyTag);
        record.id = getID(header, record.dbType);
        record.geneName = getGeneName(header, record.dbType);
        
        return record;
    }
    
    public static void main(String[] args) {
        String[] testHeaders = {
            "sp|P12345|PROT_HUMAN Protein description GN=BRCA1 PE=1 SV=2",
            "NP_001234567.1 tumor protein p53 [Homo sapiens] GN=TP53",
            "ENSP00000269305 pep chromosome:GRCh38:17:43044295:43125483:1 gene:ENSG00000141510",
            "nxp|NX_P04637|TP53 Tumor protein p53",
            "rev_sp|P12345|PROT_HUMAN Decoy protein GN=BRCA1"
        };
        
        for (String header : testHeaders) {
            ProteinRecord record = processHeader(header, "rev_");
            System.out.println("Header: " + record.originalHeader);
            System.out.println("  Type: " + record.dbType);
            System.out.println("  ID: " + record.id);
            System.out.println("  Gene: " + record.geneName);
            System.out.println("  Decoy: " + record.isDecoy);
            System.out.println();
        }
    }
}

