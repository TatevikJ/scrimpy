#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def process_input_file(ranks_bed_path, filtered_cell_barcodes_path, ranks_filtered_out_path, ratio_out_path, ratio_normalised_out_path):    
    filtered_cell_barcodes = pd.read_csv(filtered_cell_barcodes_path, 
                                         sep='\t',
                                         header=None)

    ranks_df = pd.read_csv(ranks_bed_path, 
                           sep='\t',
                           header=None)
    # Add column names
    ranks_df.columns = ['fragment_chrom', 'fragment_chromStart', 'fragment_chromEnd', 'fragment_barcode', 'fragment_readSupport',
                        'bin_chrom', 'bin_chromStart', 'bin_chromEnd', 'bin_rt_value', 'bin_rank', 'overlap_size']

    # FILTER FRAGMENT RANKS DF--------------------------------------------------------------
    # Filter fragment ranks df to take only fragments from filtered list of cell barcodes
    ranks_df = ranks_df[ranks_df["fragment_barcode"].isin(filtered_cell_barcodes[0])] 

    # Filter fragment ranks df to take only unique fragments taking max overlap size in case of duplicates
    duplicated_mask = ranks_df.duplicated(subset=["fragment_chrom", "fragment_chromStart", "fragment_chromEnd", "fragment_barcode"],
                                                  keep=False)
    
    ranks_df_duplicated = ranks_df[duplicated_mask]
    ranks_df = ranks_df[~duplicated_mask] #note: overwriting ranks_df 

    # Take first occurence of max overlap_size
    max_idx = ranks_df_duplicated.groupby(['fragment_chrom', 'fragment_chromStart', 'fragment_chromEnd', 'fragment_barcode'])['overlap_size'].idxmax()
    ranks_df_duplicated_filtered = ranks_df_duplicated.loc[max_idx] 

    # Get final df by concatenating filtered duplicates df and unique df
    ranks_df = pd.concat([ranks_df, ranks_df_duplicated_filtered]).sort_values(by=['fragment_chrom', 'fragment_chromStart']).reset_index(drop=True)

    # Save the filtered final df of fragment ranks 
    ranks_df.to_csv(ranks_filtered_out_path, 
                    sep="\t",
                    header=False,
                    index=False) 

    # GET CELL BY RANK DF------------------------------------------------------------------
    # Get rank percentages per cell and the ratio of first 3 rank percentages' sum over last 3 rank percentages' sum
    
    # Pivot ranks df to get cell (rows) by rank (columns) df
    rank_counts_per_cell = ranks_df[['fragment_barcode', 'bin_rank']].pivot_table(index='fragment_barcode', 
                                                                                  columns='bin_rank', 
                                                                                  aggfunc='size', 
                                                                                  fill_value=0)
    # Sum up each row's all values to make a new column 'total'
    rank_counts_per_cell['total'] = rank_counts_per_cell.sum(axis=1)
    # Divide each row's values by that row's total value times 100 to get % of each rank in a given cell
    rank_counts_per_cell = rank_counts_per_cell.div(rank_counts_per_cell['total'], axis=0)*100
    # Remove 'total' column
    rank_counts_per_cell = rank_counts_per_cell.drop(columns=['total'])
    # Sum up each rows 1,2,3 columns and divide by sum of that row's 8,9,10 colums to get  early/late ratio column
    rank_counts_per_cell['early_late_ratio'] = rank_counts_per_cell[[1,2,3]].sum(axis=1)/(rank_counts_per_cell[[8,9,10]].sum(axis=1))

    # Save cell by ranks df
    rank_counts_per_cell.reset_index().to_csv(ratio_out_path, 
                                              sep="\t",
                                              index=False) 

    # NORMALISE EARLY/LATE RATIO-------------------------------------------------------------
    # Normalise rank percentages of all cells by mean of rank percentages of 8% cells having lowest early/late ratio
    
    # Get mean of rank percentages of those cells having early_late_ratio <= 8% percentile value representing G1/M/G2 cells 
    # NOTE: removing nas for percentile calculation
    mean_percentage_rank_percentile_group = rank_counts_per_cell[rank_counts_per_cell.early_late_ratio <= np.percentile(rank_counts_per_cell.early_late_ratio.replace([np.inf, -np.inf], np.nan).dropna(), 8)].mean()
    # Normalise all rank percentages of all cells by obtained mean of rank percentages 
    rank_counts_per_cell_transformed = rank_counts_per_cell.div(mean_percentage_rank_percentile_group, axis=1)
    # Sum up each rows 1,2,3 normalised columns and divide by sum of that row's 8,9,10 normalised colums to get new early_late_ratio column
    rank_counts_per_cell_transformed['early_late_ratio'] = rank_counts_per_cell_transformed[[1,2,3]].sum(axis=1)/(rank_counts_per_cell_transformed[[8,9,10]].sum(axis=1))

    # Save normalised cell by rank df
    rank_counts_per_cell_transformed.reset_index().to_csv(ratio_normalised_out_path, 
                                                              sep="\t",
                                                              index=False)

if __name__ == "__main__":
    # Read input arguments from command-line
    ranks_bed_path = sys.argv[1]
    filtered_cell_barcodes_path = sys.argv[2]
    ranks_filtered_out_path = sys.argv[3]
    ratio_out_path = sys.argv[4]
    ratio_normalised_out_path = sys.argv[5]

    # Process the input files
    process_input_file(ranks_bed_path, filtered_cell_barcodes_path, ranks_filtered_out_path, ratio_out_path, ratio_normalised_out_path)