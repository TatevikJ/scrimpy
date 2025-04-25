#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def process_input_file(rt_tsv_path, out_path):    
    # 3 columns: chromosome, bin start position, rt value
    rt_df = pd.read_csv(rt_tsv_path, sep='\t')
    # Rename columns
    rt_df.columns = ['bin_chrom', 'bin_chromStart', 'bin_rt_value']

    # Bin the population RT values with percentiles and label from 10 to 1
    percentiles = np.percentile(rt_df.iloc[:,2], [i for i in range(0, 110, 10)])
    rt_df['bin_rank'] = pd.cut(rt_df.iloc[:,2], 
                              bins=percentiles, 
                              labels=[i for i in range(10, 0, -1)], 
                              include_lowest=True)

    

    # Add column for bin end position
    rt_df['bin_chromEnd'] = rt_df.iloc[:,1] + 50000

    # Save resulting df 
    rt_df.filter(['bin_chrom', 'bin_chromStart', 'bin_chromEnd', 'bin_rt_value', 'bin_rank']).to_csv(out_path,
                                                                                                  sep="\t",
                                                                                                  header=False,
                                                                                                  index=False) 


if __name__ == "__main__":
    # Read input arguments from command-line
    rt_tsv_path = sys.argv[1]
    out_path = sys.argv[2]
    
    # Process the input file
    process_input_file(rt_tsv_path, out_path)