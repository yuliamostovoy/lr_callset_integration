import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Symbolic records without SVLEN in INFO are kept unchanged.
 */
public class AddSvlenToSymbolicAlt {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i;
        int nRecords, nEdited, svlen;
        String str, alt, info, svlenStr;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; nEdited=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            alt=tokens[4]; info=tokens[7];
            if (alt.charAt(0)!='<') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            svlenStr=getInfoField(info,"SVLEN");
            if (svlenStr!=null) {
                svlen=Integer.parseInt(svlenStr);
                if (info.indexOf("SVTYPE=DEL")>=0) { 
                    alt="<DEL-"+svlen+">";
                    nEdited++;
                }
                else if (info.indexOf("SVTYPE=INV")>=0) { 
                    alt="<INV-"+svlen+">";
                    nEdited++;
                }
                else if (info.indexOf("SVTYPE=DUP")>=0) { 
                    alt="<DUP-"+svlen+">";
                    nEdited++;
                }
                else if (info.indexOf("SVTYPE=INS")>=0) { 
                    alt="<INS-"+svlen+">";
                    nEdited++;
                }
                else if (info.indexOf("SVTYPE=CNV")>=0) {
                    alt="<CNV-"+svlen+">";
                    nEdited++;
                }
                tokens[4]=alt;
            }
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
        System.err.println("nEdited="+nEdited+" ("+((100.0*nEdited)/nRecords)+"% of nRecords)");
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