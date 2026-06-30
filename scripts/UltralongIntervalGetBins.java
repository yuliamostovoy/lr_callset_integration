import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only interval calls (e.g. DEL, INV, DUP, but not 
 * INS or BND), the program prints an output BED (zero-based, left-inclusive,
 * right-exclusive) with N+4 bins for each record: N equal-sized partitions of
 * the interval, one bin before and one bin after the interval, and two bins
 * centered at the interval breakpoints.
 */
public class UltralongIntervalGetBins {
    
    private static HashMap<String,Integer> fai;
    
    /**
     * Output format: CHROM,START,END,VCFID_BINID
     *
     * @param args
     * 2: if zero, the program prints only the two bins centered at the 
     *    breakpoints;
     * 3: fixed length of each bin centered at a breakpoint.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String INPUT_FAI = args[1];
        final int N_BINS = Integer.parseInt(args[2]);
        final int BREAKPOINT_BIN_LENGTH = Integer.parseInt(args[3]);
        
        int i;
        int chromLength;
        long p, pos, svlen, quantum, binStart, binEnd;  // Long needed to avoid overflow
        String str, chrom, id, info;
        BufferedReader br;
        String[] tokens;
        
        loadFai(INPUT_FAI);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); quantum=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            chrom=tokens[0];
            chromLength=fai.get(chrom).intValue();
            pos=Integer.parseInt(tokens[1]);  // 1-based, exclusive.
            id=tokens[2];
            info=tokens[7];
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            if (N_BINS>0) {
                quantum=svlen/N_BINS;
                p=pos-quantum;
                System.out.println(chrom+"\t"+(p>=0?p:0)+"\t"+(p+quantum<=chromLength?p+quantum:chromLength)+"\t"+id+"_before");
            }
            p=pos-BREAKPOINT_BIN_LENGTH/2;
            System.out.println(chrom+"\t"+(p>=0?p:0)+"\t"+(p+BREAKPOINT_BIN_LENGTH<=chromLength?p+BREAKPOINT_BIN_LENGTH:chromLength)+"\t"+id+"_left");
            if (N_BINS>0) {
                for (i=0; i<N_BINS; i++) {
                    binStart=pos+(i*svlen)/N_BINS; binEnd=pos+((i+1)*svlen)/N_BINS;
                    System.out.println(chrom+"\t"+(binStart>=0?binStart:0)+"\t"+(binEnd<=chromLength?binEnd:chromLength)+"\t"+id+"_bin"+i);
                }
            }
            p=pos+svlen-BREAKPOINT_BIN_LENGTH/2;
            System.out.println(chrom+"\t"+(p>=0?p:0)+"\t"+(p+BREAKPOINT_BIN_LENGTH<=chromLength?p+BREAKPOINT_BIN_LENGTH:chromLength)+"\t"+id+"_right");
            if (N_BINS>0) {
                p=pos+svlen;
                System.out.println(chrom+"\t"+(p>=0?p:0)+"\t"+(p+quantum<=chromLength?p+quantum:chromLength)+"\t"+id+"_after");
            }
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
    }
    
    
    private static final void loadFai(String path) throws IOException {
        String str;
        BufferedReader br;
        String[] tokens;
        
        fai = new HashMap<String,Integer>();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            fai.put(tokens[0],Integer.valueOf(tokens[1]));
            str=br.readLine();
        }
        br.close();
    }
    
    
	/**
	 * @return NULL if $field$ does not occur in $info$.
	 */
	private static final String getInfoField(String info, String field) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return null;
		q=info.indexOf(";",p+FIELD_LENGTH);
		return info.substring(p+FIELD_LENGTH,q<0?info.length():q);
	}

}