import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF with only ultralong INS and a VCF with only ultralong DUP, the 
 * program finds every INS that has a corresponding DUP (i.e. there is a DUP 
 * that contains the INS and that has similar length).
 * 
 * Remark: length similarity is computed as in `truvari bench`.
 * 
 * Remark: this is quadratic just for simplicity, could be made much faster.
 */
public class UltralongMatchInsToDup {
    
    /**
     * @param args
     * 5: only INS of at most this length are considered.
     */
    public static void main(String[] args) throws IOException {
        final String INS_VCF_GZ = args[0];
        final String DUP_VCF_GZ = args[1];
        final int DUP_VCF_NRECORDS = Integer.parseInt(args[2]);
        final double PCTSIZE = Double.parseDouble(args[3]);
        final int SLACK_BP = Integer.parseInt(args[4]);
        final int MAX_SVLEN = Integer.parseInt(args[5]);
;
        boolean found;
        int i;
        int nIns, nMatches;
        double pos, svlen;
        String str, chr, info, value;
        BufferedReader br;
        double[][] dupPosLen;
        String[] dupChr, tokens;

        // Loading the DUP file in memory
        dupPosLen = new double[DUP_VCF_NRECORDS][3];
        dupChr = new String[DUP_VCF_NRECORDS];
        br = new BufferedReader( new InputStreamReader( (DUP_VCF_GZ.length()>=7&&DUP_VCF_GZ.substring(DUP_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(DUP_VCF_GZ)) : new FileInputStream(DUP_VCF_GZ) ) );
        str=br.readLine(); i=-1;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            i++;
            tokens=str.split("\t");
            chr=tokens[0];
            pos=Double.parseDouble(tokens[1]);  // 1-based, exclusive.
            info=tokens[7];
            svlen=Double.parseDouble(getInfoField(info,"SVLEN"));
            dupChr[i]=chr; dupPosLen[i][0]=pos; dupPosLen[i][1]=svlen;
            value=getInfoField(info,"INS_ALT");
            if (value!=null) dupPosLen[i][2]=value.length();
            else dupPosLen[i][2]=dupPosLen[i][1];
            str=br.readLine();
        }
        br.close();

        // Filtering the INS file
        br = new BufferedReader( new InputStreamReader( (INS_VCF_GZ.length()>=7&&INS_VCF_GZ.substring(INS_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INS_VCF_GZ)) : new FileInputStream(INS_VCF_GZ) ) );
        str=br.readLine(); nIns=0; nMatches=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nIns++;
            tokens=str.split("\t");
            chr=tokens[0];
            pos=Double.parseDouble(tokens[1]);  // 1-based, exclusive.
            info=tokens[7];
            svlen=Double.parseDouble(getInfoField(info,"SVLEN"));
            if (svlen>MAX_SVLEN) {
                str=br.readLine();
                continue;
            }
            found=false;
            for (i=0; i<DUP_VCF_NRECORDS; i++) {
                if (dupChr[i].equals(chr) && pos>=dupPosLen[i][0]-SLACK_BP && pos<=dupPosLen[i][0]+dupPosLen[i][1]+SLACK_BP && Math.min(svlen,dupPosLen[i][2])/Math.max(svlen,dupPosLen[i][2])>=PCTSIZE) {
                    found=true;
                    break;
                }
            }
            if (found) {
                nMatches++;
                System.out.print(tokens[0]);
                for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
                System.out.println();
            }
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println(nMatches+" INS out of "+nIns+" match a DUP ("+((100.0*nMatches)/nIns)+"%)");
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