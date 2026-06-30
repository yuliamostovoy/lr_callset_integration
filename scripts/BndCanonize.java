import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;


/**
 * In AoU, a single-sample BND VCF that results from the truvari collapse of the
 * 3 callers contains only Sniffles and pbsv calls.
 * 
 * - Pbsv emits two BND records per breakend, with correct MATEID fields.
 * - Sniffles emits one BND record per breakend (obviously without MATEID); such
 *   a record is not always the one from, say, the smallest chromosome, e.g.:
 * 
 *   chr20   29256203     id      N       N]chr4:188558260]
 * 
 * - Intra-sample truvari collapse is run in such a way that pbsv records are
 *   preferred over Sniffles records.
 *  
 * It follows that, if truvari collapse merges a pbsv and a Sniffles record,
 * the result of the merge has a paired record in output even though the merge
 * happened on just one side; but if a Sniffles record is not merged with any
 * pbsv record, the result is asymmetric. Thus, if we then perform inter-sample 
 * `bcftools merge` and `truvari collapse`, the result remains asymmetric.
 * Moreover, inter-sample collapse may not merge the same event across samples 
 * just because different samples represent it from different sides: this is 
 * because both `bcftools merge` and `truvari collapse` work on local windows of
 * the ref. (of course we could manually post-process the output, but it is
 * inelegant).
 * 
 * SVIM-asm outputs two BND records per breakend (without MATEID). An asymmetric
 * Sniffles record might match a symmetric SVIM-asm record: this is ok, but it 
 * introduces an imbalance between TPs from pbsv (represented twice) and TPs 
 * from Sniffles (represented once). Of course this is a general problem that
 * affects all records, not just TPs.
 * 
 * This program makes sure that every breakend is represented by just one BND
 * record in canonical form. This should solve all issues above and might give a
 * speedup in annotation/scoring.
 * 
 * Remarks: 
 * 1. the output VCF is not necessarily sorted;
 * 2. for simplicity, a symmetrized record uses N in REF and ALT;
 * 3. for simplicity, BNDs that do not follow the simple form (without inserted
 *    sequence) are discarded; no such BND occurs in the 292 HPRC Y2 samples at
 *    15x;
 * 4. SVIM-asm's truth should also be canonized.
 */
public class BndCanonize {
    /**
     * Not using `_` since it can appear in contig names.
     */
    private static final char CANONICAL_SEPARATOR = '@';
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        int nPaired, nUnpaired, error, nRecordsInput, nRecordsOutput, nErrorUnrecognized, nErrorSame, nErrorInsertion;
        String str, key, value;
        StringBuilder sb;
        HashMap<String,String> canonized;
        BufferedReader br;
        String[] tokens;
        
        // 1. Outputting paired records in canonical form
        canonized = new HashMap<String,String>();
        sb = new StringBuilder();
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        nRecordsInput=0; nRecordsOutput=0; nPaired=0; nUnpaired=0; nErrorUnrecognized=0; nErrorSame=0; nErrorInsertion=0;
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) System.out.println("##INFO=<ID=SYMMETRIZED,Number=0,Type=Flag,Description=\"The BND record has been symmetrized to represent it in canonical form\">");
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecordsInput++;
            tokens=str.split("\t");
            error=canonize(tokens[0],Integer.parseInt(tokens[1]),tokens[4],sb);
            if (error==1) nErrorUnrecognized++;
            else if (error==2) nErrorInsertion++;
            else if (error==3) nErrorSame++;
            else {
                key=sb.toString(); value=canonized.get(key);
                if (value==null) canonized.put(key,str);
                else {
                    nPaired+=2;
                    if (isCanonical(tokens,key)) System.out.println(str);
                    else System.out.println(value);
                    nRecordsOutput++;
                    canonized.remove(key);
                }
            }
            str=br.readLine();
        }
        br.close();
        nUnpaired=canonized.size();
        if (nRecordsInput!=nPaired+nUnpaired+nErrorUnrecognized+nErrorSame+nErrorInsertion) {
            System.err.println("ERROR: nRecordsInput = "+nRecordsInput+" != nPaired+nUnpaired+nErrorUnrecognized+nErrorSame+nErrorInsertion = "+(nPaired+nUnpaired+nErrorUnrecognized+nErrorSame+nErrorInsertion));
            System.exit(1);
        }

        // 2. Outputting unpaired records in canonical form
        canonized.forEach((k,v) -> {
            String[] tokensPrime = v.split("\t");
            if (isCanonical(tokensPrime,k)) System.out.println(v);
            else { symmetrize(tokensPrime); System.out.println(String.join("\t",tokensPrime)); }
        });
        nRecordsOutput+=nUnpaired;
        if (nRecordsOutput!=((nPaired/2)+nUnpaired)) {
            System.err.println("ERROR: nRecordsOutput = "+nRecordsOutput+" != nPaired/2+nUnpaired = "+((nPaired/2)+nUnpaired));
            System.exit(1);
        }

        // Basic counts
        System.err.println(nErrorUnrecognized+" records discarded because of unrecognized ALT ("+((nErrorUnrecognized*100.0)/nRecordsInput)+"%).");
        System.err.println(nErrorSame+" records discarded because of same CHROM,POS at their endpoints ("+((nErrorSame*100.0)/nRecordsInput)+"%).");
        System.err.println(nErrorInsertion+" records discarded because of BND with insertion ("+((nErrorInsertion*100.0)/nRecordsInput)+"%).");
        System.err.println(nPaired+" paired input records ("+((nPaired*100.0)/nRecordsInput)+"%).");
        System.err.println(nUnpaired+" unpaired input records ("+((nUnpaired*100.0)/nRecordsInput)+"%).");
        System.err.println(nRecordsInput+" total input records");
        System.err.println(nRecordsOutput+" total output records ("+((nRecordsOutput*100.0)/nRecordsInput)+"%).");
    }


    /**
     * Stores in `out` a representation of the input BND in canonical form.
     * 
     * @return
     * 0: success;
     * 1: unrecognized BND ALT;
     * 2: BND has insertion;
     * 3: REF and ALT have the same CHROM and POS.
     */
    private static final int canonize(String refChrom, int refPos, String alt, StringBuilder out) {
        char separator, refDirection, altDirection;
        int p, q;
        int first, altPos;
        String altChrom;

        // Extracting key quantities
        refDirection = (alt.charAt(0)!='[' && alt.charAt(0)!=']') ? '1':'0';  // 1 = Left
        altDirection = alt.indexOf(']')>=0 ? '1':'0';                         // 1 = Left
        p=alt.indexOf('['); q=alt.indexOf(']'); first=-1; separator='_';
        if (p>=0) { separator='['; first=p; }
        else if (q>=0) { separator=']'; first=q; }
        else return 1;
        if (p>1 || q>1) return 2;
        p=alt.indexOf(':',first+1);
        altChrom=alt.substring(first+1,p);
        q=alt.indexOf(separator,p+1);
        if (q<alt.length()-2) return 2;
        altPos=Integer.parseInt(alt.substring(p+1,q));
        
        // Canonizing
        out.delete(0,out.length());
        p=refChrom.compareTo(altChrom);
        if (p<0) { 
            out.append(refChrom); out.append(CANONICAL_SEPARATOR); out.append(refPos); out.append(CANONICAL_SEPARATOR); out.append(refDirection); out.append(CANONICAL_SEPARATOR); 
            out.append(altChrom); out.append(CANONICAL_SEPARATOR); out.append(altPos); out.append(CANONICAL_SEPARATOR); out.append(altDirection);
            return 0;
        }
        else if (p>0) {
            out.append(altChrom); out.append(CANONICAL_SEPARATOR); out.append(altPos); out.append(CANONICAL_SEPARATOR); out.append(altDirection); out.append(CANONICAL_SEPARATOR); 
            out.append(refChrom); out.append(CANONICAL_SEPARATOR); out.append(refPos); out.append(CANONICAL_SEPARATOR); out.append(refDirection);
            return 0;
        }
        if (refPos<altPos) {
            out.append(refChrom); out.append(CANONICAL_SEPARATOR); out.append(refPos); out.append(CANONICAL_SEPARATOR); out.append(refDirection); out.append(CANONICAL_SEPARATOR); 
            out.append(altChrom); out.append(CANONICAL_SEPARATOR); out.append(altPos); out.append(CANONICAL_SEPARATOR); out.append(altDirection);
            return 0;
        }
        else if (refPos>altPos) {
            out.append(altChrom); out.append(CANONICAL_SEPARATOR); out.append(altPos); out.append(CANONICAL_SEPARATOR); out.append(altDirection); out.append(CANONICAL_SEPARATOR); 
            out.append(refChrom); out.append(CANONICAL_SEPARATOR); out.append(refPos); out.append(CANONICAL_SEPARATOR); out.append(refDirection);
            return 0;
        }
        else return 3;
    }


    /**
     * @param key assumed to be the canonized form of `tokens`, and to have
     * different CHROM,POS at its endpoints;
     * @return TRUE iff the current record agrees with its canonical form.
     */
    private static final boolean isCanonical(String[] tokens, String key) {
        int p, q;

        p=key.indexOf(CANONICAL_SEPARATOR);
        if (key.substring(0,p).equals(tokens[0])) {
            q=key.indexOf(CANONICAL_SEPARATOR,p+1);
            return key.substring(p+1,q).equals(tokens[1]); 
        }
        else return false;
    }


    /**
     * Symmetrizes a BND record stored in `tokens` by changing only
     * CHROM,POS,REF,ALT and adding the SYMMETRIZED flag to INFO.
     */
    private static final void symmetrize(String[] tokens) {
        boolean refDirection, altDirection;
        char c, separator;
        int p, q;
        int first;
        String refChrom, refPos, alt, altChrom, altPos;

        refChrom=tokens[0]; refPos=tokens[1]; alt=tokens[4];

        // Extracting key quantities
        c=alt.charAt(0);
        refDirection=(c!='[')&&(c!=']');   // True = Left
        altDirection=alt.indexOf(']')>=0;  // True = Left
        p=alt.indexOf('['); q=alt.indexOf(']'); first=-1; separator='_';
        if (p>=0) { separator='['; first=p; }
        else if (q>=0) { separator=']'; first=q; }
        else {
            System.err.println("ERROR: unrecognized ALT = "+alt);
            System.exit(1);
        }
        p=alt.indexOf(':',first+1);
        altChrom=alt.substring(first+1,p);
        q=alt.indexOf(separator,p+1);
        altPos=alt.substring(p+1,q);

        // Symmetrizing
        tokens[0]=altChrom; tokens[1]=altPos; tokens[3]="N";
        separator=refDirection?']':'[';
        if (altDirection) tokens[4]="N"+separator+refChrom+":"+refPos+separator;
        else tokens[4]=separator+refChrom+":"+refPos+separator+"N";
        tokens[7]=tokens[7].equals(".")?"SYMMETRIZED":tokens[7]+";SYMMETRIZED";
    }

}