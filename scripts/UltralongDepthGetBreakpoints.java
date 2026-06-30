import java.util.Arrays;
import java.util.Vector;
import java.io.*;


/**
 * Given the output of `samtools depth` over a region, the program computes a
 * simple estimate of the leftmost and rightmost positions that delimit an 
 * increase in coverage.
 */
public class UltralongDepthGetBreakpoints {
    
    /**
     * @args
     */
    public static void main(String[] args) throws IOException {
        final String SAMTOOLS_DEPTH_FILE = args[0];
        final int SAMTOOLS_DEPTH_N_POSITIONS = Integer.parseInt(args[1]);
        final int BIN_LENGTH = Integer.parseInt(args[2]);
        final double BIN_COVERAGE_RATIO = Double.parseDouble(args[3]);

        final int MAX_BREAKPOINTS = 30;  // Arbitrary
        
        int i, j, p;
        int length, sum1, sum2, leftBreakpoint, rightBreakpoint, posStart, posEnd, validPairs, maxLength, bestLeft, bestRight;
        long firstPos;
        double lowDepth1, lowDepth2, ratio, maxRatio;
        String str;
        BufferedReader br;
        int[] depth;
        Vector<Breakpoint> leftBreakpoints, rightBreakpoints;

        // Loading all depths in memory
        depth = new int[SAMTOOLS_DEPTH_N_POSITIONS]; firstPos=-1;
        br = new BufferedReader(new FileReader(SAMTOOLS_DEPTH_FILE));
        str=br.readLine(); i=-1;
        while (str!=null) {
            length=str.length();
            p=0;
            while (p<length) {
                if (str.charAt(p)=='\t') break;
                else p++;
            }
            p++; posStart=p;
            while (p<length) {
                if (str.charAt(p)=='\t') break;
                else p++;
            }
            posEnd=p;
            if (firstPos==-1) firstPos=Integer.parseInt(str.substring(posStart,posEnd));
            p++;
            depth[++i]=Integer.parseInt(str.substring(p));
            str=br.readLine();
        }
        br.close();

        // Collecting all left breakpoints, if any.
        leftBreakpoints = new Vector<Breakpoint>();
        sum1=0;
        for (i=0; i<BIN_LENGTH; i++) sum1+=depth[i];
        sum2=0;
        for (i=BIN_LENGTH; i<2*BIN_LENGTH; i++) sum2+=depth[i];
        if (sum1>=BIN_LENGTH) {
            maxRatio=((double)sum2)/sum1;
            if (maxRatio>=BIN_COVERAGE_RATIO) leftBreakpoints.add(new Breakpoint(i-BIN_LENGTH,((double)sum1)/BIN_LENGTH,maxRatio));   
        }
        else maxRatio=0;
        for (i=2*BIN_LENGTH; i<SAMTOOLS_DEPTH_N_POSITIONS; i++) {
            sum1+=depth[i-BIN_LENGTH]-depth[i-2*BIN_LENGTH];
            sum2+=depth[i]-depth[i-BIN_LENGTH];
            if (sum1<BIN_LENGTH) continue;
            ratio=((double)sum2)/sum1;
            if (ratio>maxRatio) maxRatio=ratio;
            if (ratio>=BIN_COVERAGE_RATIO) leftBreakpoints.add(new Breakpoint(i-BIN_LENGTH,((double)sum1)/BIN_LENGTH,ratio));
        }
        if (leftBreakpoints.isEmpty()) return;
        System.err.println(leftBreakpoints.size()+" left breakpoints, maxRatio="+String.format("%.2f",maxRatio));

        // Collecting all right breakpoints, if any.
        rightBreakpoints = new Vector<Breakpoint>();
        sum2=0;
        for (i=SAMTOOLS_DEPTH_N_POSITIONS-1; i>=SAMTOOLS_DEPTH_N_POSITIONS-BIN_LENGTH; i--) sum2+=depth[i];
        sum1=0;
        for (i=SAMTOOLS_DEPTH_N_POSITIONS-BIN_LENGTH-1; i>=SAMTOOLS_DEPTH_N_POSITIONS-2*BIN_LENGTH; i--) sum1+=depth[i];
        if (sum2>=BIN_LENGTH) {
            maxRatio=((double)sum1)/sum2;
            if (maxRatio>=BIN_COVERAGE_RATIO) rightBreakpoints.add(new Breakpoint(i+BIN_LENGTH,((double)sum2)/BIN_LENGTH,maxRatio));
        }
        else maxRatio=0;
        for (i=SAMTOOLS_DEPTH_N_POSITIONS-2*BIN_LENGTH-1; i>=0; i--) {
            sum1+=depth[i]-depth[i+BIN_LENGTH];
            sum2+=depth[i+BIN_LENGTH]-depth[i+2*BIN_LENGTH];
            if (sum2<BIN_LENGTH) continue;
            ratio=((double)sum1)/sum2;
            if (ratio>maxRatio) maxRatio=ratio;
            if (ratio>=BIN_COVERAGE_RATIO) rightBreakpoints.add(new Breakpoint(i+BIN_LENGTH,((double)sum2)/BIN_LENGTH,ratio));
        }
        if (rightBreakpoints.isEmpty()) return;
        System.err.println(rightBreakpoints.size()+" right breakpoints, maxRatio="+String.format("%.2f",maxRatio));

        // Picking the top-k breakpoints if there are too many
        if (leftBreakpoints.size()>MAX_BREAKPOINTS) {
            leftBreakpoints.sort(null);
            leftBreakpoints.setSize(MAX_BREAKPOINTS);
        }
        if (rightBreakpoints.size()>MAX_BREAKPOINTS) {
            rightBreakpoints.sort(null);
            rightBreakpoints.setSize(MAX_BREAKPOINTS);
        }
        
        // Finding a longest valid pair, if any.
        // Remark: medians are computed naively and should be made faster.
        validPairs=0; maxLength=0; bestLeft=-1; bestRight=-1; maxRatio=0;
        for (i=0; i<leftBreakpoints.size(); i++) {
            leftBreakpoint=leftBreakpoints.get(i).breakpoint;
            lowDepth1=leftBreakpoints.get(i).lowDepth;
            for (j=0; j<rightBreakpoints.size(); j++) {
                rightBreakpoint=rightBreakpoints.get(j).breakpoint;
                if (rightBreakpoint<=leftBreakpoint) continue;
                lowDepth2=rightBreakpoints.get(j).lowDepth;
                int[] newArray = Arrays.copyOfRange(depth,leftBreakpoint,rightBreakpoint+1);
                Arrays.sort(newArray);
                if (newArray[newArray.length/2]>=(lowDepth1<lowDepth2?lowDepth1:lowDepth2)*BIN_COVERAGE_RATIO) {
                    validPairs++;
                    if (rightBreakpoint-leftBreakpoint>maxLength) {
                        maxLength=rightBreakpoint-leftBreakpoint;
                        bestLeft=leftBreakpoint; bestRight=rightBreakpoint;
                    }
                }
                ratio=((double)newArray[newArray.length/2])/(lowDepth1<lowDepth2?lowDepth1:lowDepth2);
                if (ratio>maxRatio) maxRatio=ratio;
            }
        }
        System.err.println("validPairs="+validPairs+" maxLength="+maxLength+" maxRatio="+String.format("%.2f",maxRatio));
        if (validPairs==0) return;

        // Outputting
        System.out.println((firstPos+bestLeft)+"\t"+(firstPos+bestRight));
    }


    private static class Breakpoint implements Comparable {
        int breakpoint;
        double lowDepth, ratio;

        public Breakpoint(int breakpoint, double lowDepth, double ratio) {
            this.breakpoint=breakpoint;
            this.lowDepth=lowDepth;
            this.ratio=ratio;
        }

        public boolean equals(Object other) {
            Breakpoint otherBreakpoint = (Breakpoint)other;
            return this.ratio==otherBreakpoint.ratio;
        }

        /**
         * Decreasing ratio
         */
        public int compareTo(Object other) {
            Breakpoint otherBreakpoint = (Breakpoint)other;
            if (this.ratio>otherBreakpoint.ratio) return -1;
            else if (this.ratio<otherBreakpoint.ratio) return 1;
            else return 0;
        }
    }

}