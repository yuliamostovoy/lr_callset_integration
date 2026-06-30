import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only svim-asm INS records, the program separates
 * records that are likely DUP based on the dipcall confident BED.
 * 
 * Remark: the confident BED is assumed to contain sorted, non-overlapping
 * intervals.
 * 
 * Remark: the output DUP VCF is not necessarily sorted if it contains DUP 
 * records.
 */
public class UltralongSvimasmInsExtractDups {
    
    /**
     * @param args
     * 7: 0=output INS->DUP records as DUP; 1=output INS->DUP records as INS.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String INPUT_DIPCALL_BED = args[1];
        final int N_BED_RECORDS = Integer.parseInt(args[2]);
        final int SLACK_BP = Integer.parseInt(args[3]);
        final double LENGTH_SIMILARITY = Double.parseDouble(args[4]);
        final String OUTPUT_VCF_DUP = args[5];
        final String OUTPUT_VCF_INS = args[6];
        final int OUTPUT_DUP_MODE = Integer.parseInt(args[7]);

        boolean found;
        int i, j;
        int pos, newPos, newEnd, newLength, bedLength, svlen, nRecords, nDups;
        String str, chrom, alt, info;
        BufferedReader br;
        BufferedWriter bwIns, bwDup;
        int[] chrFirst;
        String[] bedChr, tokens;
        int[][] bedIntervals;

        // Loading the entire BED in memory
        bedChr = new String[N_BED_RECORDS];
        bedIntervals = new int[N_BED_RECORDS][2];
        i=-1;
        br = new BufferedReader(new FileReader(INPUT_DIPCALL_BED));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            i++;
            bedChr[i]=tokens[0];
            bedIntervals[i][0]=Integer.parseInt(tokens[1]);
            bedIntervals[i][1]=Integer.parseInt(tokens[2]);
            str=br.readLine();
        }
        br.close();
        chrFirst = new int[25];
        Arrays.fill(chrFirst,-1);
        for (i=0; i<N_BED_RECORDS; i++) {
            j=chrom2index(bedChr[i]);
            if (j!=-1 && chrFirst[j]==-1) chrFirst[j]=i;
        }

        // Processing the VCF
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
            chrom=tokens[0];
            pos=Integer.parseInt(tokens[1]);
            info=tokens[7];
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            i=chrom2index(chrom);
            if (i==-1) { 
                // Non-standard chromosome
                bwIns.write(str+"\n");
                str=br.readLine();
                continue; 
            }
            i=chrFirst[i]; 
            if (i==-1) { 
                // No BED record for this chromosome
                bwIns.write(str+"\n");
                str=br.readLine();
                continue; 
            }
            found=false;
            while (i<N_BED_RECORDS) {
                if (!bedChr[i].equals(chrom)) break;
                else if (pos<bedIntervals[i][0]-SLACK_BP) break;
                else if (pos>bedIntervals[i][1]+SLACK_BP) { i++; continue; }
                else if (svlen>=(bedIntervals[i][1]-bedIntervals[i][0])*LENGTH_SIMILARITY) {
                    found=true; nDups++;
                    if (OUTPUT_DUP_MODE==0) {
                        newPos=bedIntervals[i][0]; newEnd=bedIntervals[i][1]; newLength=newEnd-newPos;
                        tokens[1]=newPos+"";
                        alt=tokens[4]; tokens[4]="<DUP>";
                        info=addOrReplaceInfoField(info,"SVLEN",String.valueOf(newLength));
                        info=addOrReplaceInfoField(info,"SVTYPE","DUP");
                        info=addOrReplaceInfoField(info,"END",newEnd+"");
                        info+=";INS_POS="+pos+";INS_ALT="+alt+";INS_QUAL="+tokens[5];
                        tokens[7]=info;
                    }
                    bwDup.write(tokens[0]);
                    for (j=1; j<tokens.length; j++) bwDup.write("\t"+tokens[j]);
                    bwDup.write("\n");
                    break;
                }
                else { i++; continue; }
            }
            if (!found) bwIns.write(str+"\n");
            str=br.readLine();
        }
        br.close(); bwDup.close(); bwIns.close();
        System.err.println("Extracted "+nDups+" DUPs out of "+nRecords+" INS records ("+String.format("%.2f",(100.0*nDups)/nRecords)+"%) with dipcall BED.");
    }


    /**
     * @return 0-based. -1 if `chrom` is not a standard chromosome.
     */
    private static final int chrom2index(String chrom) {
        final char c = chrom.charAt(3);
        if (c=='X' || c=='x') return 22;
        else if (c=='Y' || c=='y') return 23;
        else if (c=='M' || c=='m') return 24;
        else { try { return Integer.parseInt(chrom.substring(3))-1; } catch (NumberFormatException e) { return -1; } }
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