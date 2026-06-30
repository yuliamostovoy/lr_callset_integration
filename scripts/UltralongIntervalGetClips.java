import java.util.*;
import java.io.*;


/**
 * Given a SAM file and a window, the program outputs two files containing a
 * record for every alignment that has a left (resp. right) soft- or hard-clip
 * at a position inside the window. We call such alignments left- (resp. right-)
 * maximal.
 *
 * Output format: READ_ID,IS_RC,READ_POS,READ_LENGTH
 */
public class UltralongIntervalGetClips {
    
    /**
     * @param args 
     * 1: 0-based, inclusive;
     * 2: 0-based, exclusive.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        final int START = Integer.parseInt(args[1]);
        final int END = Integer.parseInt(args[2]);
        final int MIN_CLIP_LENGTH = Integer.parseInt(args[3]);
        final String OUTPUT_PREFIX = args[4];
        
        final int RC_MASK = 16;
        
        boolean isRc, matchFound;
        char c;
        int i, j;
        int refPos, readPos, newReadPos, cigarLength, readLength, length;
        String str, cigar, output;
        BufferedReader br;
        BufferedWriter bwLeft, bwRight;
        String[] tokens;
        
        bwLeft = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_leftmaximal.txt"));
        bwRight = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_rightmaximal.txt"));
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            refPos=Integer.parseInt(tokens[3]);  // 1-based, first mapping pos.
            cigar=tokens[5];
            cigarLength=cigar.length();
            isRc=(Integer.parseInt(tokens[1])&RC_MASK)!=0;
            readLength=tokens[9].length();
            readPos=0;  // 0-based, first mapping pos.
            i=0; matchFound=false;
            for (j=1; j<cigarLength; j++) {
                c=cigar.charAt(j);
                if (c=='M' || c=='=' || c=='X') {
                    matchFound=true;
                    length=Integer.parseInt(cigar.substring(i,j));
                    refPos+=length; readPos+=length;
                    i=j+1;
                }
                else if (c=='I') {
                    readPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='D') {
                    matchFound=true;
                    refPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='N') {
                    refPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='S') {
                    newReadPos=readPos+Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                    if ( refPos-1>=START && refPos-1<END &&
                         ( (matchFound && readLength-readPos>=MIN_CLIP_LENGTH) || (!matchFound && newReadPos>=MIN_CLIP_LENGTH) )
                       ) {
                        output=tokens[0]+"\t"+(isRc?"1\t":"0\t")+(matchFound?readPos:newReadPos)+"\t"+(readLength+"\n");
                        if (!matchFound) bwLeft.write(output);
                        else bwRight.write(output);
                    }
                    readPos=newReadPos;
                }
                else if (c=='H') {
                    newReadPos=readPos+Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                    if ( refPos-1>=START && refPos-1<END &&
                         ( (matchFound && readLength-readPos>=MIN_CLIP_LENGTH) || (!matchFound && newReadPos>=MIN_CLIP_LENGTH) )
                       ) {
                        output=tokens[0]+"\t"+(isRc?"1\t":"0\t")+(matchFound?readPos:newReadPos)+"\t"+(readLength+"\n");
                        if (!matchFound) bwLeft.write(output);
                        else bwRight.write(output);
                    }
                }
                else if (c=='P') i=j+1;
            }

            // Next iteration
            str=br.readLine();
        }
        br.close(); bwLeft.close(); bwRight.close();
    }

}