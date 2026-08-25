# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
import pysam
import os
import json
from tqdm import tqdm
from scipy.stats import mannwhitneyu, ks_2samp
import warnings
warnings.filterwarnings('ignore')
import sys




FEATURES = [
    'log_svlen', 'depth_ratio', 'depth_mad', 'ab', 'cn_slop',
    'mq_drop', 'clip_frac', 'split_reads', 'read_len_med', 'strand_bias',
    'gc_frac', 'homopolymer_max', 'lcr_mask',
    'support_read', 'svtype_DEL'
]
     

def check_existing_data(dataset_name):
    """Check for existing processed data"""
    complete_file = f'./complete_features.csv'

    if os.path.exists(complete_file):
        try:
            existing_data = pd.read_csv(complete_file)
            existing_variants = set()
            for _, row in existing_data.iterrows():
                variant_id = f"{row.chrom}:{row.pos}:{row.end}:{row.label}"
                existing_variants.add(variant_id)

            print(f"Found existing data for {dataset_name}: {len(existing_data):,} variants")
            return existing_data, existing_variants
        except Exception as e:
            print(f"Error reading existing file: {e}")

    print(f"No existing data found for {dataset_name}")
    return None, set()


def check_dataset_files(dataset_name, dataset_info):
    """Check if all required files exist"""
    print(f"\nChecking {dataset_name}:")

    required_files = ['bam', 'tp_comp_vcf', 'fp_vcf', 'ref']
    all_exist = True

    for file_type in required_files:
        filepath = dataset_info[file_type]
        exists = os.path.exists(filepath)
        status = "OK" if exists else "MISSING"
        print(f"  {file_type}: {status}")
        if not exists:
            all_exist = False

    return all_exist


def load_variants_with_filtering(vcf_path, label, existing_variants):
    """Load variants, filtering out already processed ones"""
    variants = []
    vcf = VCF(vcf_path)

    total_count = 0
    new_count = 0

    for variant in vcf:
        total_count += 1

        var_info = {
            'label': label,
            'chrom': variant.CHROM,
            'pos': variant.POS,
            'end': variant.end if hasattr(variant, 'end') else variant.POS + 1,
            'svlen': variant.INFO.get('SVLEN', 0),
            'svtype': variant.INFO.get('SVTYPE', 'UNK'),
            'id': variant.ID
        }

        variant_id = f"{var_info['chrom']}:{var_info['pos']}:{var_info['end']}:{label}"

        if variant_id not in existing_variants:
            variants.append((variant, var_info))
            new_count += 1

    vcf.close()
    print(f"    {label}: {new_count:,} new / {total_count:,} total")
    return variants
     

def compute_size_features(bam_file, ref_file, chrom, start, end, svlen, svtype):
    """Compute size and copy number features"""
    features = {}

    # log_svlen
    abs_svlen = abs(svlen) if svlen != 0 else 1
    features['log_svlen'] = np.log10(abs_svlen)

    try:
        # Handle insertions differently
        if svtype == 'INS':
            window_size = max(100, abs(svlen) // 10)
            sv_start = start - window_size // 2
            sv_end = start + window_size // 2
            sv_length = window_size
        else:
            sv_start = start
            sv_end = end
            sv_length = end - start

        flank = 1000

        # Get reads
        sv_reads = list(bam_file.fetch(chrom, sv_start, sv_end))
        left_reads = list(bam_file.fetch(chrom, max(0, sv_start - flank), sv_start))
        right_reads = list(bam_file.fetch(chrom, sv_end, sv_end + flank))

        # Calculate coverages
        sv_coverage = len(sv_reads) / sv_length if sv_length > 0 else 0
        left_coverage = len(left_reads) / flank
        right_coverage = len(right_reads) / flank
        control_coverage = (left_coverage + right_coverage) / 2

        # depth_ratio
        features['depth_ratio'] = sv_coverage / control_coverage if control_coverage > 0 else np.nan

        # depth_mad
        coverage_array = np.zeros(sv_length)
        for read in sv_reads:
            read_start = max(0, read.reference_start - sv_start)
            read_end = min(sv_length, read.reference_end - sv_start)
            if read_end > read_start:
                coverage_array[read_start:read_end] += 1

        if len(coverage_array) > 0 and sv_length > 0:
            median_cov = np.median(coverage_array)
            features['depth_mad'] = np.median(np.abs(coverage_array - median_cov)) / sv_coverage if sv_coverage > 0 else np.nan
        else:
            features['depth_mad'] = np.nan

        # ab (allele balance)
        total_coverage = sv_coverage + control_coverage
        features['ab'] = sv_coverage / total_coverage if total_coverage > 0 else np.nan

        # cn_slop (distal comparison)
        try:
            distal_start = max(0, start - 25000)
            distal_end = start - 5000
            if distal_end > distal_start:
                distal_reads = list(bam_file.fetch(chrom, distal_start, distal_end))
                distal_coverage = len(distal_reads) / (distal_end - distal_start)
                features['cn_slop'] = sv_coverage / distal_coverage if distal_coverage > 0 else np.nan
            else:
                features['cn_slop'] = np.nan
        except:
            features['cn_slop'] = np.nan

    except Exception as e:
        features.update({
            'depth_ratio': np.nan, 'depth_mad': np.nan,
            'ab': np.nan, 'cn_slop': np.nan
        })

    return features
     

def compute_read_quality_features(bam_file, chrom, start, end, svtype):
    """Compute read quality and mapping features"""
    features = {}

    try:
        # Handle insertions differently
        if svtype == 'INS':
            window_size = 200
            sv_start = start - window_size // 2
            sv_end = start + window_size // 2
        else:
            sv_start = start
            sv_end = end

        # Get reads
        sv_reads = list(bam_file.fetch(chrom, sv_start, sv_end))
        left_reads = list(bam_file.fetch(chrom, max(0, start-1000), start))
        right_reads = list(bam_file.fetch(chrom, end, end+1000))

        # mq_drop
        sv_mq = [read.mapping_quality for read in sv_reads if read.mapping_quality is not None]
        flank_mq = [read.mapping_quality for read in left_reads + right_reads if read.mapping_quality is not None]

        sv_mq_median = np.median(sv_mq) if sv_mq else 0
        flank_mq_median = np.median(flank_mq) if flank_mq else 0
        features['mq_drop'] = flank_mq_median - sv_mq_median

        # clip_frac, split_reads, read_len_med
        clipped_reads = 0
        split_read_count = 0
        read_lengths = []

        for read in sv_reads:
            # Check clipping
            if read.cigartuples:
                for op, length in read.cigartuples:
                    if op in [4, 5] and length >= 10:
                        clipped_reads += 1
                        break

            # Check split reads
            if read.has_tag('SA') or read.is_supplementary:
                split_read_count += 1

            # Read length
            if read.query_length is not None:
                read_lengths.append(read.query_length)

        features['clip_frac'] = clipped_reads / len(sv_reads) if sv_reads else 0
        features['split_reads'] = split_read_count / len(sv_reads) if sv_reads else 0
        features['read_len_med'] = np.median(read_lengths) if read_lengths else np.nan

        # strand_bias
        forward_reads = sum(1 for read in sv_reads if not read.is_reverse)
        reverse_reads = len(sv_reads) - forward_reads

        if len(sv_reads) > 0:
            forward_frac = forward_reads / len(sv_reads)
            reverse_frac = reverse_reads / len(sv_reads)
            features['strand_bias'] = abs(forward_frac - reverse_frac)
        else:
            features['strand_bias'] = np.nan

    except Exception as e:
        features.update({
            'mq_drop': np.nan, 'clip_frac': np.nan, 'split_reads': 0,
            'read_len_med': np.nan, 'strand_bias': np.nan
        })

    return features
     

def compute_sequence_context_features(ref_file, chrom, start, end):
    """Compute sequence context features"""
    features = {}

    try:
        # Get sequence with flanks
        flank = 500
        seq_start = max(0, start - flank)
        seq_end = end + flank
        
        # Handle chromosome naming variations
        sequence = ""
        test_chroms = [
            chrom
        ]        
        for test_chrom in test_chroms:
            if test_chrom in ref_file.references:
                sequence = ref_file.fetch(test_chrom, seq_start, seq_end).upper()
                break
        
        if sequence:
            # gc_frac
            gc_count = sequence.count('G') + sequence.count('C')
            features['gc_frac'] = gc_count / len(sequence)

            # homopolymer_max
            max_homopolymer = 0
            current_base = ''
            current_count = 0

            for base in sequence:
                if base == current_base:
                    current_count += 1
                else:
                    max_homopolymer = max(max_homopolymer, current_count)
                    current_base = base
                    current_count = 1
            max_homopolymer = max(max_homopolymer, current_count)
            features['homopolymer_max'] = max_homopolymer / len(sequence)

            # lcr_mask
            distinct_bases = len(set(sequence))
            features['lcr_mask'] = 1 if distinct_bases <= 2 else 0
        else:
            features.update({
                'gc_frac': np.nan, 'homopolymer_max': np.nan, 'lcr_mask': np.nan
            })

    except Exception as e:
        features.update({
            'gc_frac': np.nan, 'homopolymer_max': np.nan, 'lcr_mask': np.nan
        })

    return features


def compute_all_features_for_variant(variant, bam_file, ref_file):
    """Compute all 15 features for a single variant"""
    # Basic variant info
    chrom = variant.CHROM
    start = variant.POS
    variantid = variant.ID
    end = variant.end if hasattr(variant, 'end') and variant.end else start + 1
    svlen = variant.INFO.get('SVLEN', end - start)
    svtype = variant.INFO.get('SVTYPE', 'UNK')

    # Initialize features
    features = {
        'chrom': chrom, 'pos': start, 'end': end,
        'svlen': svlen, 'svtype': svtype,
        'id': variantid
    }

    # Compute feature categories
    features.update(compute_size_features(bam_file, ref_file, chrom, start, end, svlen, svtype))
    features.update(compute_read_quality_features(bam_file, chrom, start, end, svtype))
    features.update(compute_sequence_context_features(ref_file, chrom, start, end))

    # Caller features
    features['support_read'] = variant.INFO.get('SUPPORT', np.nan)
    features['svtype_DEL'] = 1 if svtype == 'DEL' else 0

    return features
     

def process_dataset(dataset_name, dataset_info, output_file):
    """Process a single dataset to extract features"""
    print(f"\nProcessing {dataset_name}")
    print("=" * 50)

    # Check existing data
    existing_data, existing_variants = check_existing_data(dataset_name)

    # Load new variants only
    print("Loading variants...")
    new_variants = load_variants_with_filtering(
        dataset_info['vcf_gz'], 'TP', existing_variants
    )

    if len(new_variants) == 0:
        print("All variants already processed")
        return existing_data

    print(f"Processing {len(new_variants):,} new variants...")

    # Open files
    try:
        bam_file = pysam.AlignmentFile(dataset_info['bam'], 'rb')
        ref_file = pysam.FastaFile(dataset_info['ref'])
    except Exception as e:
        print(f"Error opening files: {e}")
        return existing_data

    # Process new variants
    new_results = []
    failed_count = 0

    with tqdm(total=len(new_variants), desc="Computing features") as pbar:
        for i, (variant, var_info) in enumerate(new_variants):
            try:
                features = compute_all_features_for_variant(variant, bam_file, ref_file)
                features.update(var_info)
                features['dataset'] = dataset_name
                new_results.append(features)

                pbar.set_postfix({'Success': f"{len(new_results)/(i+1)*100:.1f}%"})
            except Exception as e:
                failed_count += 1
                pbar.set_postfix({'Failed': failed_count})

            pbar.update(1)

    # Close files
    #bam_file.close()
    #ref_file.close()

    # Combine data
    if new_results:
        new_df = pd.DataFrame(new_results)

        if existing_data is not None:
            combined_df = pd.concat([existing_data, new_df], ignore_index=True)
            print(f"Combined: {len(existing_data):,} existing + {len(new_df):,} new = {len(combined_df):,} total")
        else:
            combined_df = new_df
            print(f"New dataset: {len(new_df):,} variants")

        # Save complete dataset
        combined_df.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")

        return combined_df
    else:
        print("No new variants processed successfully")
        return existing_data     


def main(args):
    DATASETS = {
        'dataset': {
            'vcf_gz': sys.argv[1],
            'bam': sys.argv[2],
            'ref': sys.argv[3]
        }
    }
    output_file = sys.argv[4]
    
    all_dataframes = []
    for dataset_name, dataset_info in DATASETS.items():
        try:
            df = process_dataset(dataset_name, dataset_info, output_file)
            if df is not None and len(df) > 0:
                all_dataframes.append(df)
        except Exception as e:
            print(f"Error with {dataset_name}: {e}")

    if all_dataframes:
        # Combine all data
        combined_df = pd.concat(all_dataframes, ignore_index=True)

        print(f"\nFinal dataset summary:")
        print(f"Total variants: {len(combined_df):,}")
        print(f"Datasets: {list(combined_df['dataset'].unique())}")
        print(f"Variants: {len(combined_df[combined_df['label'] == 'TP']):,}")

        # Feature availability
        print(f"\nFeature availability:")
        for feature in FEATURES:
            if feature in combined_df.columns:
                available = combined_df[feature].notna().sum()
                total = len(combined_df)
                pct = (available / total) * 100
                print(f"  {feature:<18} {available:>7,}/{total:<7,} ({pct:>5.1f}%)")

        # Save final dataset
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        #main_file = f'./complete_features_{timestamp}.csv'
        #combined_df.to_csv(main_file, index=False)
        #print(f"\nSaved complete dataset: {main_file}")
    else:
        print("No data processed!")




if __name__ == "__main__":
    main(sys.argv)
