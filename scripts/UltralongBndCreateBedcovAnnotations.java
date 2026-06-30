import java.util.*;
import java.io.*;


/**
 * Given the BED file created by `samtools bedcov` (with format CHROM,START,END,
 * RECORDID,BINID,BEDCOV), the program reformats it as RECORDID,BINID,BEDCOV.
 * Output values are normalized by BIN_LENGTH.
 */
public class UltralongBndCreateBedcovAnnotations {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        final int BIN_LENGTH = Integer.parseInt(args[1]);
        
        int i, p;
        String str, id, recordId, binId;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_BED)));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            recordId=tokens[3];
            binId=tokens[4];
            System.out.printf("%s\t%s\t%.3f\n",recordId,binId,(Double.parseDouble(tokens[5])/BIN_LENGTH));
            str=br.readLine();
        }
        br.close();
    }

}