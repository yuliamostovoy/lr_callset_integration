import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only svim-asm INS records annotated with `truvari 
 * anno remap`, the program separates records that are classified as DUP,
 * assuming that `INFO/remap_coords` contains the accurate span of the entire
 * (possibly complex) duplication.
 * 
 * Remark: the output DUP VCF is not necessarily sorted if it contains DUP 
 * records.
 */
public class UltralongSvimasmInsExtractDupsPrime {
    
    /**
     * @param args
     * 3: 0=output INS->DUP records as DUP; 1=output INS->DUP records as INS.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_VCF_DUP = args[1];
        final String OUTPUT_VCF_INS = args[2];
        final int OUTPUT_DUP_MODE = Integer.parseInt(args[3]);

        int i, p, q;
        int pos, start, end, nRecords, nDups;
        String str, alt, info, classification, coords;
        BufferedReader br;
        BufferedWriter bwIns, bwDup;
        String[] tokens;

        bwDup = new BufferedWriter(new FileWriter(OUTPUT_VCF_DUP));
        bwIns = new BufferedWriter(new FileWriter(OUTPUT_VCF_INS));
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; nDups=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM") && OUTPUT_DUP_MODE==0) {
                    bwDup.write("##INFO=<ID=INS_POS,Number=1,Type=Integer,Description=\"The POS of the original INS record\">\n");
                    bwDup.write("##INFO=<ID=INS_ALT,Number=1,Type=String,Description=\"The ALT allele of the original INS record\">\n");
                    bwDup.write("##INFO=<ID=INS_QUAL,Number=1,Type=Integer,Description=\"The QUAL of the original INS record\">\n");
                }
                bwDup.write(str+"\n");
                bwIns.write(str+"\n");
                str=br.readLine();
                continue;
            }
            nRecords++;
            tokens=str.split("\t");
            info=tokens[7];
            classification=getInfoField(info,"remap_classification");
            if (classification.equals("tandem") || classification.equals("tandem_complex") || classification.equals("tandem_inverted")) {
                nDups++;
                if (OUTPUT_DUP_MODE==0) {
                    coords=getInfoField(info,"remap_coords");
                    p=coords.indexOf(":");
                    q=coords.indexOf("-");
                    start=Integer.parseInt(coords.substring(p+1,q));
                    end=Integer.parseInt(coords.substring(q+1));
                    pos=Integer.parseInt(tokens[1]); tokens[1]=start+"";
                    alt=tokens[4]; tokens[4]="<DUP>";
                    info=addOrReplaceInfoField(info,"END",end+"");
                    info=addOrReplaceInfoField(info,"SVLEN",(end-start)+"");
                    info=addOrReplaceInfoField(info,"SVTYPE","DUP");
                    info+=";INS_POS="+pos+";INS_ALT="+alt+";INS_QUAL="+tokens[5];
                    tokens[7]=info;
                }
                bwDup.write(tokens[0]);
                for (i=1; i<tokens.length; i++) bwDup.write("\t"+tokens[i]);
                bwDup.write("\n");
            }
            else bwIns.write(str+"\n");
            str=br.readLine();
        }
        br.close(); bwDup.close(); bwIns.close();
        System.err.println("Extracted "+nDups+" DUPs out of "+nRecords+" INS records ("+String.format("%.2f",(100.0*nDups)/nRecords)+"%) with truvari anno remap.");
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


    private static final String addOrReplaceInfoField(String info, String field, String newValue) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        if (info.equals(".")) return field+"="+newValue;
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return info+";"+field+"="+newValue;
		q=info.indexOf(";",p+FIELD_LENGTH);
        return info.substring(0,p+FIELD_LENGTH)+newValue+(q>=0?info.substring(q):"");
	}

}