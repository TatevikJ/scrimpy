#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def process_input_file(overlapping_fragments_bed_path, filtered_cell_barcodes_path, atac_fragments_tsv_gz_path, ranks_filtered_path, ratio_path, ratio_normalised_path, metrics_out_path):    
    filtered_cell_barcodes = pd.read_csv(filtered_cell_barcodes_path, 
                                         sep='\t',
                                         header=None)
    overlapping_fragments_df = pd.read_csv(overlapping_fragments_bed_path, 
                                        sep='\t',
                                        header=None)

    # Add column names
    overlapping_fragments_df.columns = ['fragment_chrom', 'fragment_chromStart', 'fragment_chromEnd', 
                                        'fragment_barcode', 'fragment_readSupport']

    # FILTER OVERLAPPING FRAGMENTS DF--------------------------------------------------------------
    # Filter overlapping fragments df to take only those from filtered list of cell barcodes
    overlapping_fragments_df_filtered = overlapping_fragments_df[overlapping_fragments_df.fragment_barcode.isin(filtered_cell_barcodes[0])]

    # GET PER CELL SUMMARY METRICS TABLE-----------------------------------------------------------
    # Read and filter atac_fragments tsv file by list of barcodes
    # atac_fragments = pd.read_csv(atac_fragments_tsv_gz_path,
    #                              compression='gzip',
    #                              sep = '\t',
    #                              header=None,
    #                              comment='#') #to ignore first commented rows
    if atac_fragments_tsv_gz_path.endswith('.gz'):
        atac_fragments = pd.read_csv(atac_fragments_tsv_gz_path,
                                      compression='gzip',
                                      sep='\t',
                                      header=None,
                                      comment='#')
    else:
        atac_fragments = pd.read_csv(atac_fragments_tsv_gz_path,
                                      sep='\t',
                                      header=None,
                                      comment='#')
    
    atac_fragments.columns = ['chrom', 'chromStart', 'chromEnd', 'barcode', 'readSupport']
    atac_fragments_filtered = atac_fragments[atac_fragments.barcode.isin(filtered_cell_barcodes[0])]

    # Read filtered fragment ranks' df
    ranks_df = pd.read_csv(ranks_filtered_path,
                          sep = '\t',
                          header = None)

    ranks_df.columns = ['fragment_chrom', 'fragment_chromStart', 'fragment_chromEnd', 'fragment_barcode', 'fragment_readSupport',
                        'bin_chrom', 'bin_chromStart', 'bin_chromEnd', 'bin_rt_value', 'bin_rank', 'overlap_size']


    # Read non-normalised and normalised rank counts and ratio per cell dfs
    rank_counts_per_cell = pd.read_csv(ratio_path, 
                                       sep="\t") 
    rank_counts_per_cell_transformed = pd.read_csv(ratio_normalised_path, 
                                                   sep="\t") 

    # Create a summary table with all metrics per cell barcode
    metrics_per_barcode = ranks_df.groupby('fragment_barcode').size().reset_index(name='num_of_background_fragments') 
    
    metrics_per_barcode = metrics_per_barcode.merge(ranks_df.groupby('fragment_barcode')['bin_rank'].mean().reset_index(name='mean_rank'), 
                                                    on = 'fragment_barcode')
    metrics_per_barcode = metrics_per_barcode.merge(ranks_df.groupby('fragment_barcode')['bin_rank'].median().reset_index(name='median_rank'), 
                                                    on = 'fragment_barcode')
    
    metrics_per_barcode = metrics_per_barcode.merge(atac_fragments_filtered.groupby('barcode').size().reset_index(name='num_of_all_fragments'), 
                                                    left_on = 'fragment_barcode', 
                                                    right_on = 'barcode', 
                                                    how='outer')
    # COUNT NUMBER OF OVERLAPPING FRAGMENTS PER CELL
    metrics_per_barcode = metrics_per_barcode.merge(overlapping_fragments_df_filtered.groupby('fragment_barcode').size().reset_index(name='num_of_overlapping_fragments'), 
                                                    left_on = 'fragment_barcode', 
                                                    right_on = 'fragment_barcode', 
                                                    how='outer')

    
    # Replace nan with 0 in num_of_overlapping_fragments column
    metrics_per_barcode['num_of_overlapping_fragments'] = metrics_per_barcode['num_of_overlapping_fragments'].replace(np.nan, 0)

    # NORMALISE OVERLAPPING FRAGMENTS COUNTS BY NUMBER OF FRAGMENTS PER CELL
    metrics_per_barcode['percentage_of_overlapping_fragments'] = metrics_per_barcode.num_of_overlapping_fragments/metrics_per_barcode.num_of_all_fragments*100

    # Add early/late ratio
    metrics_per_barcode = metrics_per_barcode.merge(rank_counts_per_cell[['fragment_barcode', 'early_late_ratio']], 
                                                    on = 'fragment_barcode', 
                                                    how = 'outer')
    metrics_per_barcode.rename(columns={"early_late_ratio": "early_late_ratio_original"}, 
                               inplace=True)
    
    # Add normalised early/late ratio
    metrics_per_barcode = metrics_per_barcode.merge(rank_counts_per_cell_transformed[['fragment_barcode', 'early_late_ratio']], 
                                                    on = 'fragment_barcode', 
                                                    how = 'outer')
    metrics_per_barcode.rename(columns={"early_late_ratio": "early_late_ratio_transformed"}, 
                               inplace=True)

    metrics_per_barcode.rename(columns={"early_late_ratio_original": "early_late_ratio"}, 
                               inplace=True)

    # metrics_per_barcode = metrics_per_barcode[metrics_per_barcode['percentage_of_overlapping_fragments'] != 0]
    
    # ASSIGN CELL CYCLE STATE BASED ON THRESHOLD
    assignment_conditions = [((metrics_per_barcode['early_late_ratio_transformed'] < 1.5) & (np.log2(metrics_per_barcode['percentage_of_overlapping_fragments']) < -1)),
                             (np.log2(metrics_per_barcode['percentage_of_overlapping_fragments']) > 1)]
    states = ['G1', 'G2M']
    metrics_per_barcode['state'] = np.select(assignment_conditions, states, default='S')

    # ASSIGN CELL CYCLE STATE USING SEURAT?
    # ...
    
    # Save metrics df
    metrics_per_barcode.to_csv(metrics_out_path,
                               index=False)
    
if __name__ == "__main__":
    # Read input arguments from command-line
    overlapping_fragments_bed_path = sys.argv[1]
    filtered_cell_barcodes_path = sys.argv[2]
    atac_fragments_tsv_gz_path = sys.argv[3]
    ranks_filtered_path = sys.argv[4]
    ratio_path = sys.argv[5]
    ratio_normalised_path = sys.argv[6]
    metrics_out_path = sys.argv[7]

    # Process the input files
    process_input_file(overlapping_fragments_bed_path, filtered_cell_barcodes_path, atac_fragments_tsv_gz_path, ranks_filtered_path, ratio_path, ratio_normalised_path, metrics_out_path)