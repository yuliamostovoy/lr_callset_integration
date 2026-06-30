import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * The program converts every INSDUP record back to INS, saving the breakpoints 
 * of the INSDUP interval in INFO.
 * 
 * Remark: the output VCF is not necessarily sorted.
 */
public class UltralongInsdups2Ins {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        int i;
        int pos, insPos, insLen, svlen, insQual, nRecords, nInsdup;
        String str, insAlt, info, insLen_str;
        BufferedReader br;
        String[] tokens;

        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; nInsdup=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) {
                    System.out.println("##INFO=<ID=INSDUP_POS,Number=1,Type=Integer,Description=\"The POS value of the INSDUP interval\">");
                    System.out.println("##INFO=<ID=INSDUP_SVLEN,Number=1,Type=Integer,Description=\"The SVLEN value of the INSDUP interval\">");
                }
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            tokens=str.split("\t");
            info=tokens[7];
            if (info.indexOf("INSDUP")<0) System.out.println(str);
            else {
                nInsdup++;
                pos=Integer.parseInt(tokens[1]);
                svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
                insPos=Integer.parseInt(getInfoField(info,"INS_POS"));
                insAlt=getInfoField(info,"INS_ALT");
                insLen_str=getInfoField(info,"INS_SVLEN");
                if (insLen_str!=null) insLen=Integer.parseInt(insLen_str);
                else insLen=insAlt.length();
                insQual=Integer.parseInt(getInfoField(info,"INS_QUAL"));
                tokens[1]=insPos+"";
                tokens[4]=insAlt;
                tokens[5]=insQual+"";
                info=deleteInfoField(info,"INS_POS");
                info=deleteInfoField(info,"INS_ALT");
                info=deleteInfoField(info,"INS_QUAL");
                if (insLen_str!=null) info=deleteInfoField(info,"INS_SVLEN");
                info=addOrReplaceInfoField(info,"SVTYPE","INS");
                info=addOrReplaceInfoField(info,"SVLEN",insLen+"");
                info=addOrReplaceInfoField(info,"END",insPos+"");
                info+=";INSDUP_POS="+pos+";INSDUP_SVLEN="+svlen;
                tokens[7]=info;
                System.out.print(tokens[0]);
                for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
                System.out.println();
            }
            str=br.readLine();
        }
        br.close();
        System.err.println(nInsdup+" INSDUP->INS conversions out of "+nRecords+" total records ("+String.format("%.2f",(100.0*nInsdup)/nRecords)+"%).");
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


    private static final String deleteInfoField(String info, String field) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        if (info.equals(".")) return info;
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return info;
		q=info.indexOf(";",p+FIELD_LENGTH);
        return info.substring(0,p)+(q>=0?info.substring(q+1):"");
	}

}