import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Similar to `CleanRefAltQual.java`.
 */
public class CleanRefAltQualBnd {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String FORCE_QUAL = args[1];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i;
        int nRecords;
        String str;
        StringBuilder buffer;
        BufferedReader br;
        String[] tokens;
        
        buffer = new StringBuilder();
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            tokens[3]=tokens[3].toUpperCase();
            tokens[4]=capitalizeBnd(tokens[4]);
            tokens[5]=FORCE_QUAL;
            tokens[6]="PASS";
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
    }


    /**
     * Capitalizes just the DNA characters
     */
    private static final String capitalizeBnd(String alt) {
        char separator;
        int p, q, first;

        p=alt.indexOf('['); q=alt.indexOf(']'); first=-1; separator='_';
        if (p>=0) { separator='['; first=p; }
        else if (q>=0) { separator=']'; first=q; }
        else {
            System.err.println("ERROR: unrecognized BND ALT: "+alt);
            System.exit(1);
        }
        p=alt.indexOf(separator,first+1);
        return alt.substring(0,first).toUpperCase()+alt.substring(first,p+1)+alt.substring(p+1).toUpperCase();
    }
    
}