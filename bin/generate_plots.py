#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def generate_plots(rt_tsv_path, metrics_csv_path, figures_out_dir):
    
    # FIGURE 1--------------------------------------------------------------
    # Note: will make same plot for each sample
    # Read replication timing rank df
    rt_df = pd.read_csv(rt_tsv_path, 
                        sep='\t',
                        header=None)
    
    percentiles = np.percentile(rt_df.iloc[:,3], [i for i in range(0, 110, 10)])

    # Plot
    p1 = plt.figure()
    plt.hist(rt_df.iloc[:,3], bins=300)
    plt.xlabel('Replication timing value')
    plt.ylabel('Frequency')
    # plt.title('Replication rank groups defined by percentiles of replication timing values')
    for p in percentiles:
        plt.axvline(x=p, color='red', linestyle='--', linewidth=1)

    # Save pdf 
    p1.savefig(figures_out_dir + '/distribution_of_replication_timing_values.pdf')


    # FIGURE 2--------------------------------------------------------------
    # Read df
    metrics_per_barcode = pd.read_csv(metrics_csv_path)

    # Plot
    p2 = plt.figure()
    plt.hist(metrics_per_barcode['num_of_background_fragments'], 
             alpha = 0.5, 
             bins=300)
    # plt.title('Distribution of the number of background fragments per cell')
    plt.xlabel('Number of background fragments per cell')
    plt.ylabel('Frequency')
    
    # Save pdf 
    p2.savefig(figures_out_dir + '/distribution_of_number_of_background_fragments_per_cell.pdf')

    # FIGURE 3--------------------------------------------------------------
    # Plot
    p3 = plt.figure()
    plt.hist(metrics_per_barcode['mean_rank'], 
             alpha = 0.5, 
             bins=300)
    # plt.title('Distribution of the mean replication timing rank of background fragments per cell')
    plt.xlabel('Mean rank per cell')
    plt.ylabel('Frequency')
    
    # Save pdf 
    p3.savefig(figures_out_dir + '/distribution_of_mean_replication_timing_rank_of_background_fragments_per_cell.pdf')

    # FIGURE 4--------------------------------------------------------------
    # Plot
    # Note: rows with early/late ratio = inf not plotted
    p4 = plt.figure()
    plt.hist(metrics_per_barcode['early_late_ratio'].replace([np.inf, -np.inf], np.nan), 
             alpha = 0.5, 
             bins=300)
    plt.axvline(x=1, color='red', linestyle='--', linewidth=0.7)
    plt.text(x=1, y=-30, s=f'{1}', color='r', horizontalalignment='center', fontsize=8)
    plt.axvline(x=2, color='red', linestyle='--', linewidth=0.7)
    plt.text(x=2, y=-30, s=f'{2}', color='r', horizontalalignment='center', fontsize=8)
    # plt.title('Distribution of early/late ratio')
    plt.xlabel('Early/late ratio (not normalised)')
    plt.ylabel('Frequency')
    
    # Save pdf 
    p4.savefig(figures_out_dir + '/distribution_of_early_late_ratio.pdf')

    # FIGURE 5--------------------------------------------------------------
    # Plot
    # Note: rows with normalised early/late ratio = inf not plotted
    p5 = plt.figure()
    plt.hist(metrics_per_barcode['early_late_ratio_transformed'].replace([np.inf, -np.inf], np.nan), 
             alpha = 0.5, 
             bins=300)
    plt.axvline(x=1, color='red', linestyle='--', linewidth=0.7)
    plt.text(x=1, y=-30, s=f'{1}', color='r', horizontalalignment='center', fontsize=8)
    plt.axvline(x=2, color='red', linestyle='--', linewidth=0.7)
    plt.text(x=2, y=-30, s=f'{2}', color='r', horizontalalignment='center', fontsize=8)
    # plt.title('Distribution of early/late ratio after normalisation')
    plt.xlabel('Early/late ratio (normalised)')
    plt.ylabel('Frequency')
    
    # Save pdf 
    p5.savefig(figures_out_dir + '/distribution_of_early_late_ratio_after_normalisation.pdf')


    # # FIGURE 6--------------------------------------------------------------
    # # Plot
    # # Note: rows with 0% overlapping fragments and 0 fragment count with late rank (leading to ealry/late ratio getting inf value) not plotted
    # p6 = sns.jointplot(x=metrics_per_barcode.early_late_ratio_transformed.replace([np.inf, -np.inf], np.nan), 
    #               y=np.log2(metrics_per_barcode.percentage_of_overlapping_fragments.replace(0, np.nan)), 
    #               kind='hex')
    
    # # Adjust the figure to add a full-height color bar
    # fig = p6.fig
    # fig.subplots_adjust(right=0.85)  # Make space on the right for the color bar
    
    # # Add color bar to a new axis that spans the full height of the figure
    # cax = fig.add_axes([0.87, 0.09, 0.03, 0.76])  # Position: [left, bottom, width, height]
    # cbar = fig.colorbar(p6.ax_joint.collections[0], cax=cax, orientation='vertical')
    # cbar.set_label('Number of points')
    
    # # Set axis labels and title
    # p6.set_axis_labels('Early/late ratio (normalised)', 'Log2 percentage of overlapping fragments')
    # # p6.ax_joint.set_title('2D histogram of early/late ratio and log2 percentage of overlapping fragments')
    
    # # Add lines and texts
    # p6.ax_joint.vlines(x=1.5, ymin=-5, ymax=-1, color='red', linestyle='--', linewidth=0.7)
    # p6.ax_joint.text(x=1.5, y=-6, s=f'{1.5}', color='r', horizontalalignment='center', fontsize=6)
    
    # p6.ax_joint.hlines(y=-1, xmin=0, xmax=1.5, color='red', linestyle='--', linewidth=0.7)
    # p6.ax_joint.text(x=-0.8, y=-1, s=f'{-1}', color='r', horizontalalignment='center', fontsize=6)
    
    # p6.ax_joint.hlines(y=1, xmin=0, xmax=10, color='red', linestyle='--', linewidth=0.7)
    # p6.ax_joint.text(x=-0.8, y=1, s=f'{1}', color='r', horizontalalignment='center', fontsize=6)
    
    # p6.ax_joint.set_xticks(np.arange(0, 11))

    # # Save pdf 
    # p6.savefig(figures_out_dir + '/2D_histogram_of_early_late_ratio_and_log2_percentage_of_overlapping_fragments.pdf')



    #gpt
    # FIGURE 6--------------------------------------------------------------
    # Plot
    # Note: rows with 0% overlapping fragments and 0 fragment count with late rank (leading to early/late ratio getting inf value) not plotted
    p6 = sns.jointplot(x=metrics_per_barcode.early_late_ratio_transformed.replace([np.inf, -np.inf], np.nan), 
                       y=np.log2(metrics_per_barcode.percentage_of_overlapping_fragments.replace(0, np.nan)), 
                       kind='hex')
    
    # Adjust the figure to add a full-height color bar
    fig = p6.fig
    fig.subplots_adjust(right=0.85)  # Make space on the right for the color bar
    
    # Add color bar to a new axis that spans the full height of the figure
    cax = fig.add_axes([0.87, 0.09, 0.03, 0.76])  # Position: [left, bottom, width, height]
    cbar = fig.colorbar(p6.ax_joint.collections[0], cax=cax, orientation='vertical')
    cbar.set_label('Number of points')
    
    # Set axis labels and title
    p6.set_axis_labels('Early/late ratio (normalised)', 'Log2 percentage of overlapping fragments')
    # p6.ax_joint.set_title('2D histogram of early/late ratio and log2 percentage of overlapping fragments')

    # Get axis limits dynamically
    x_limits = p6.ax_joint.get_xlim()
    y_limits = p6.ax_joint.get_ylim()
    
    # Add vertical and horizontal lines with dynamic limits
    p6.ax_joint.vlines(x=1.5, ymin=min(-1, y_limits[0]), ymax=-1, color='red', linestyle='--', linewidth=0.7)
    p6.ax_joint.text(x=1.5, y=y_limits[0] - 0.1 * (y_limits[1] - y_limits[0]), 
                     s=f'{1.5}', color='r', horizontalalignment='center', fontsize=6)
    
    p6.ax_joint.hlines(y=-1, xmin=x_limits[0], xmax=1.5, color='red', linestyle='--', linewidth=0.7)
    p6.ax_joint.text(x=x_limits[0] - 0.1 * (x_limits[1] - x_limits[0]), y=-1, 
                     s=f'{-1}', color='r', horizontalalignment='center', fontsize=6)
    
    p6.ax_joint.hlines(y=1, xmin=x_limits[0], xmax=x_limits[1], color='red', linestyle='--', linewidth=0.7)
    p6.ax_joint.text(x=x_limits[0] - 0.1 * (x_limits[1] - x_limits[0]), y=1, 
                     s=f'{1}', color='r', horizontalalignment='center', fontsize=6)
    
    # Set xticks dynamically if needed, or leave as is
    #p6.ax_joint.set_xticks(np.arange(np.floor(x_limits[0]), np.ceil(x_limits[1]) + 1))
    
    # Save pdf 
    p6.savefig(figures_out_dir + '/2D_histogram_of_early_late_ratio_and_log2_percentage_of_overlapping_fragments.pdf')



if __name__ == "__main__":
    # Read input arguments from command-line
    rt_tsv_path = sys.argv[1]
    metrics_csv_path = sys.argv[2]
    figures_out_dir = sys.argv[3]

    # Generate plots using input arguments
    generate_plots(rt_tsv_path, metrics_csv_path,figures_out_dir)

