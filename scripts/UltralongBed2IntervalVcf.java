import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * The program transforms an input BED into a minimal headerless VCF with only
 * symbolic calls.
 */
public class UltralongBed2IntervalVcf {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        final String SVTYPE = args[1];
        
        int i;
        int start, end;
        String str, chrom;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_BED)));
        str=br.readLine(); i=-1;
        while (str!=null) {
            tokens=str.split("\t");
            chrom=tokens[0];
            start=Integer.parseInt(tokens[1]);  // 0-based, inclusive.
            end=Integer.parseInt(tokens[2]);    // 0-based, exclusive.
            System.out.println(chrom+"\t"+(start+1)+"\t"+(++i)+"\tN\t<"+SVTYPE+">\t60\tPASS\tSVTYPE="+SVTYPE+";SVLEN="+(end-start)+";END="+end+"\tGT\t0/1");
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
    }

}