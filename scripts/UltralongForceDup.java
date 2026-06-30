import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only DUP records, possibly marked with different 
 * DUP subtypes, the program forces every ALT and SVTYPE to be simply DUP.
 */
public class UltralongForceDup {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            tokens[4]="<DUP>";
            tokens[7]=addOrReplaceInfoField(tokens[7],"SVTYPE","DUP");
            System.out.print(tokens[0]);
            for (int i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
            System.out.println();
            str=br.readLine();
        }
        br.close();
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