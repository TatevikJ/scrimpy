#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Required parameters
params.input_paths_csv = false
params.gene_regions_bed = false
params.rt_tsv = false
params.outdir = false

// Optional parameters with defaults
params.mem_per_cpu = '40G'
params.cpus = 4
params.queue = 'long'

// Print parameter summary
log.info """
    SCRIMPy: Single Cell Replication Inference from Multiome data using Python
    ==========================================================================
    input_paths_csv : ${params.input_paths_csv}
    gene_regions_bed: ${params.gene_regions_bed}
    rt_tsv         : ${params.rt_tsv}
    outdir         : ${params.outdir}
    """

// Validate required inputs
if (!params.input_paths_csv || !params.gene_regions_bed || !params.rt_tsv || !params.outdir) {
    error "Missing required parameters. Please provide input_paths_csv, gene_regions_bed, rt_tsv, and outdir"
}

// Main workflow
workflow {
    sample_info_ch = Channel
                        .fromPath(params.input_paths_csv)
                        .splitCsv(header: true)
                        .map { row ->
                            // Extract the values from each row
                            def sample = row.sample
                            def peaks = row.peaks
                            def fragments = row.fragments  
                            def barcodes = row.barcodes
                            
                            // Return a tuple with the values
                            tuple(sample, file(peaks), file(fragments), file(barcodes))
                        }    

    filtered_fragments = FILTER_FRAGMENTS(sample_info_ch)

    gene_regions_ch = Channel.fromPath(params.gene_regions_bed)

    filtered_fragments_combined_ch = filtered_fragments.combine(gene_regions_ch)

    background_fragments = GET_BACKGROUND_FRAGMENTS(filtered_fragments_combined_ch)

    rt_tsv_ch = Channel.fromPath(params.rt_tsv)

    rt_ranks = GET_REPLICATION_RANKS(rt_tsv_ch)

    background_fragments_combined_ch = background_fragments.combine(rt_ranks)

    fragment_ranks = GET_FRAGMENT_RANKS(background_fragments_combined_ch)

    //processed_fragments = PROCESS_FRAGMENT_RANKS(fragment_ranks)

    overlapping_fragments = GET_OVERLAPPING_FRAGMENTS(background_fragments)

    // Join by sample name
    //overlapping_fragments_processed_fragments_joined_ch = overlapping_fragments.join(processed_fragments, by: 0)

    //overlapping_fragments_processed_fragments_joined_ch.view()

    //metrics = PROCESS_OVERLAPPING_FRAGMENTS(overlapping_fragments_processed_fragments_joined_ch)

    //metrics_combined_ch = metrics.combine(rt_ranks)

    //GENERATE_PLOTS(metrics_combined_ch)
}


// Process 1: Filter atac fragments by provided list of barcodes and keep only fragments from autosomal chromosomes for each sample
process FILTER_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(peaks), path(fragments), path(barcodes)
    
    output:
    tuple val(sample), path(peaks), path("atac_fragments_filtered.tsv.gz"), path(barcodes)
    
    script:
    """
    if [[ ${barcodes} == *.gz ]]; then
        zcat ${barcodes} | awk 'NR==FNR {filter[\$1]; next} \$4 in filter' - <(zcat ${fragments}) \\
            | awk -F'\\t' '\$1 ~ /^(chr([0-9]{1,2})|([0-9]{1,2}))\$/ {print}' \\
            > atac_fragments_filtered.tsv
    else
        zcat ${fragments} | awk 'NR==FNR {filter[\$1]; next} \$4 in filter' ${barcodes} - \\
            | awk -F'\\t' '\$1 ~ /^(chr([0-9]{1,2})|([0-9]{1,2}))\$/ {print}' \\
            > atac_fragments_filtered.tsv
    fi
    
    gzip atac_fragments_filtered.tsv
    """
}

// Process 2: Get background (exclude peaks or peaks+genes) and peak fragments for each sample
process GET_BACKGROUND_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(peaks), path(filtered_fragments), path(barcodes), path(gene_regions)
    
    output:
    tuple val(sample), 
          path("atac_background_fragments.bed"),
          path("atac_background_fragments_exclude_gene_regions.bed"),
          path("atac_peak_fragments.bed"),
          path(peaks),
          path(filtered_fragments),
          path(barcodes)
    
    script:
    """
    # Exclude peaks
    bedtools intersect -v \\
        -a ${filtered_fragments} \\
        -b ${peaks} \\
        > atac_background_fragments.bed
    
    # Exclude peaks + intragenic regions
    bedtools intersect -v \\
        -a atac_background_fragments.bed \\
        -b ${gene_regions} \\
        > atac_background_fragments_exclude_gene_regions.bed
    
    # Get peaks
    bedtools intersect -wa \\
        -a ${filtered_fragments} \\
        -b ${peaks} \\
        > atac_peak_fragments.bed
    """
}

// Process 3: Get replication ranks per bin for each sample
process GET_REPLICATION_RANKS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path rt_tsv
    
    output:
    path "rt_ranks.tsv"
    
    script:
    """
    get_replication_ranks_per_bin.py ${rt_tsv} rt_ranks.tsv
    """
}

// Process 4: Get fragment ranks for each sample
process GET_FRAGMENT_RANKS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), 
          path(background_fragments), 
          path(background_fragments_exclude_gene),
          path(peak_fragments),
          path(peaks),
          path(filtered_fragments),
          path(barcodes),
          path(rt_ranks)
    
    output:
    tuple val(sample),
          path("atac_fragments_ranks.bed"),
          path("atac_background_fragments_ranks.bed"),
          path("atac_background_fragments_exclude_gene_regions_ranks.bed"),
          path(background_fragments),
          path(background_fragments_exclude_gene),
          path(peak_fragments),
          path(barcodes)
    
    script:
    """
    # Get ranks for all fragments
    bedtools intersect -wo \\
        -a ${filtered_fragments} \\
        -b ${rt_ranks} \\
        > atac_fragments_ranks.bed
        
    # Get ranks for background fragments
    bedtools intersect -wo \\
        -a ${background_fragments} \\
        -b ${rt_ranks} \\
        > atac_background_fragments_ranks.bed
        
    # Get ranks for background fragments excluding gene regions
    bedtools intersect -wo \\
        -a ${background_fragments_exclude_gene} \\
        -b ${rt_ranks} \\
        > atac_background_fragments_exclude_gene_regions_ranks.bed
    """
}

// Process 5: Process fragment ranks: 1. filter fragment ranks, 2. get rank counts and ratio per cell, normalise ratio
// !!! Note: revise this part, include ratio and normalisation?
process PROCESS_FRAGMENT_RANKS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample),
          path(fragments_ranks),
          path(background_fragments_ranks),
          path(background_fragments_exclude_gene_regions_ranks),          
          path(background_fragments),
          path(background_fragments_exclude_gene),
          path(peak_fragments),
          path(barcodes)
    
    output:
    tuple val(sample),
          path("atac_fragments_ranks_filtered.tsv"),
          path("atac_fragments_ranks_counts_and_ratio_per_cell.tsv"),
          path("atac_fragments_ranks_counts_and_ratio_per_cell_normalised.tsv"),
          path("atac_background_fragments_ranks_filtered.tsv"),
          path("atac_background_fragments_ranks_counts_and_ratio_per_cell.tsv"),
          path("atac_background_fragments_ranks_counts_and_ratio_per_cell_normalised.tsv"),
          path("atac_background_fragments_exclude_gene_regions_ranks_filtered.tsv"),
          path("atac_background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell.tsv"),
          path("atac_background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell_normalised.tsv")

    script:
    """
    echo "CPUs: \$task.cpus"
    echo "Memory: \$task.memory"

    # Process all fragments
    process_background_fragments_ranks.py \\
        ${fragments_ranks} \\
        ${barcodes} \\
        atac_fragments_ranks_filtered.tsv \\
        atac_fragments_ranks_counts_and_ratio_per_cell.tsv \\
        atac_fragments_ranks_counts_and_ratio_per_cell_normalised.tsv

    # Process background fragments
    process_background_fragments_ranks.py \\
        ${background_fragments_ranks} \\
        ${barcodes} \\
        atac_background_fragments_ranks_filtered.tsv \\
        atac_background_fragments_ranks_counts_and_ratio_per_cell.tsv \\
        atac_background_fragments_ranks_counts_and_ratio_per_cell_normalised.tsv

    # Process background fragments excluding gene regions 
    process_background_fragments_ranks.py \\
        ${background_fragments_exclude_gene_regions_ranks} \\
        ${barcodes} \\
        atac_background_fragments_exclude_gene_regions_ranks_filtered.tsv \\
        atac_background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell.tsv \\
        atac_background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell_normalised.tsv
    """
}


// Process 6: Get overlapping fragments with different region exclusion criteria for each sample
process GET_OVERLAPPING_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), 
          path(background_fragments), 
          path(background_fragments_exclude_gene),
          path(peak_fragments),
          path(peaks),
          path(filtered_fragments),
          path(barcodes)
    
    output:
    tuple val(sample),
          path("atac_overlapping_fragments.bed"),
          path("atac_overlapping_fragments_exclude_peak_regions.bed"),
          path("atac_overlapping_fragments_exclude_peak_gene_regions.bed"),
          path("atac_overlapping_fragments_only_peak_regions.bed"),
          path(filtered_fragments),
          path(barcodes),
          path(background_fragments), 
          path(background_fragments_exclude_gene),
          path(peak_fragments)
          
    script:
    """    
    #### Get overlapping fragments - exclude nothing
    bedtools intersect \\
        -a ${filtered_fragments} \\
        -b ${filtered_fragments} \\
        -wo -sorted | \\
    awk '\$1 == \$6 && \$4 == \$9 && (\$2 < \$7 || (\$2 == \$7 && \$3 < \$8))' - \\
        > atac_overlapping_fragment_pairs.bed
    
    awk '{print \$1, \$2, \$3, \$4, \$5; print \$6, \$7, \$8, \$9, \$10}' OFS="\\t" atac_overlapping_fragment_pairs.bed | \\
    sort -k1,1 -k2,2n -k3,3n -u \\
        > atac_overlapping_fragments.bed

    #### Get overlapping fragments - exclude peaks
    bedtools intersect \\
        -a ${background_fragments} \\
        -b ${background_fragments} \\
        -wo -sorted | \\
    awk '\$1 == \$6 && \$4 == \$9 && (\$2 < \$7 || (\$2 == \$7 && \$3 < \$8))' - \\
        > atac_overlapping_fragment_pairs_exclude_peak_regions.bed
    
    awk '{print \$1, \$2, \$3, \$4, \$5; print \$6, \$7, \$8, \$9, \$10}' OFS="\\t" atac_overlapping_fragment_pairs_exclude_peak_regions.bed | \\
    sort -k1,1 -k2,2n -k3,3n -u \\
        > atac_overlapping_fragments_exclude_peak_regions.bed

    #### Get overlapping fragments - exclude peaks+intragenic regions
    bedtools intersect \
         -a ${background_fragments_exclude_gene} \\
         -b ${background_fragments_exclude_gene} \\
         -wo  \\
         -sorted | \\
	awk '\$1 == \$6 && \$4 == \$9 && (\$2 < \$7 || (\$2 == \$7 && \$3 < \$8))' - \\
	    > atac_overlapping_fragment_pairs_exclude_peak_gene_regions.bed
    
    awk '{print \$1, \$2, \$3, \$4, \$5; print \$6, \$7, \$8, \$9, \$10}' OFS="\\t" atac_overlapping_fragment_pairs_exclude_peak_gene_regions.bed | \\
	sort -k1,1 -k2,2n -k3,3n -u \\
	     > atac_overlapping_fragments_exclude_peak_gene_regions.bed

    #### Get overlapping fragments - only peak regions
    bedtools intersect \
         -a ${peak_fragments} \\
         -b ${peak_fragments} \\
         -wo  \\
         -sorted | \\
	awk '\$1 == \$6 && \$4 == \$9 && (\$2 < \$7 || (\$2 == \$7 && \$3 < \$8))' - \\
	    > atac_overlapping_fragment_pairs_only_peak_regions.bed
    
    awk '{print \$1, \$2, \$3, \$4, \$5; print \$6, \$7, \$8, \$9, \$10}' OFS="\\t" atac_overlapping_fragment_pairs_only_peak_regions.bed | \\
	sort -k1,1 -k2,2n -k3,3n -u \\
	     > atac_overlapping_fragments_only_peak_regions.bed
    """
}

// Process 7: Process overlapping fragments and generate metrics
process PROCESS_OVERLAPPING_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample),
          path(overlapping_fragments),
          path(overlapping_fragments_exclude_peak_regions),
          path(overlapping_fragments_exclude_peak_gene_regions),
          path(overlapping_fragments_only_peak_regions),
          path(filtered_fragments),
          path(barcodes),
          path(background_fragments), 
          path(background_fragments_exclude_gene),
          path(peak_fragments),
          path(fragments_ranks_filtered),
          path(fragments_ranks_counts_and_ratio_per_cell),
          path(fragments_ranks_counts_and_ratio_per_cell_normalised),
          path(background_fragments_ranks_filtered),
          path(background_fragments_ranks_counts_and_ratio_per_cell),
          path(background_fragments_ranks_counts_and_ratio_per_cell_normalised),
          path(background_fragments_exclude_gene_regions_ranks_filtered),
          path(background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell),
          path(background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell_normalised)

    output:
    tuple val(sample), 
          path("metrics_per_cell_exclude_peak_regions_only_background.csv"),
          path("metrics_per_cell_exclude_peak_gene_regions_only_background.csv"),
          path("metrics_per_cell_exclude_nothing.csv"),
          path("metrics_per_cell_exclude_peak_regions_only_background.csv"),
          path("metrics_per_cell_exclude_peak_gene_regions.csv"),
          path("metrics_per_cell_only_peak_regions.csv")

    script:
    """    
    #### Exclude peaks for early/late ratio and exclude nothing for overlapping fragments
    process_overlapping_fragments.py \\
        ${overlapping_fragments} \\
        ${barcodes} \\
        ${filtered_fragments} \\
        ${background_fragments_ranks_filtered} \\
        ${background_fragments_ranks_counts_and_ratio_per_cell} \\
        ${background_fragments_ranks_counts_and_ratio_per_cell_normalised} \\
        metrics_per_cell_exclude_peak_regions_only_background.csv

    #### Exclude peaks and gene regions for early/late ratio and exclude nothing for overlapping fragments
    process_overlapping_fragments.py \\
        ${overlapping_fragments} \\
        ${barcodes} \\
        ${filtered_fragments} \\
        ${background_fragments_exclude_gene_regions_ranks_filtered} \\
        ${background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell} \\
        ${background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell_normalised} \\
        metrics_per_cell_exclude_peak_gene_regions_only_background.csv

    #### Exclude nothing
    process_overlapping_fragments.py \\
        ${overlapping_fragments} \\
        ${barcodes} \\
        ${filtered_fragments} \\
        ${fragments_ranks_filtered} \\
        ${fragments_ranks_counts_and_ratio_per_cell} \\
        ${fragments_ranks_counts_and_ratio_per_cell_normalised} \\
        metrics_per_cell_exclude_nothing.csv

    #### Exclude peaks
    process_overlapping_fragments.py \\
        ${overlapping_fragments_exclude_peak_regions} \\
        ${barcodes} \\
        ${background_fragments} \\
        ${background_fragments_ranks_filtered} \\
        ${background_fragments_ranks_counts_and_ratio_per_cell} \\
        ${background_fragments_ranks_counts_and_ratio_per_cell_normalised} \\
        metrics_per_cell_exclude_peak_regions_only_background.csv

    #### Exclude peaks+intragenic regions
    process_overlapping_fragments.py \\
        ${overlapping_fragments_exclude_peak_gene_regions} \\
        ${barcodes} \\
        ${background_fragments_exclude_gene_regions} \\
        ${background_fragments_exclude_gene_regions_ranks_filtered} \\
        ${background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell} \\
        ${background_fragments_exclude_gene_regions_ranks_counts_and_ratio_per_cell_normalised} \\
        metrics_per_cell_exclude_peak_gene_regions.csv

    #### Exclude peaks for early/late ratio and use only peaks for overlapping fragments
    process_overlapping_fragments.py \\
        ${overlapping_fragments_only_peak_regions} \\
        ${barcodes} \\
        ${peak_fragments} \\
        ${background_fragments_ranks_filtered} \\
        ${background_fragments_ranks_counts_and_ratio_per_cell} \\
        ${background_fragments_ranks_counts_and_ratio_per_cell_normalised} \\
        metrics_per_cell_only_peak_regions.csv
    """
}

// Process 8: Generate plots
process GENERATE_PLOTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/figures", mode: 'copy'
    
    input:
    tuple val(sample), 
      path("metrics_per_cell_exclude_peak_regions_only_background.csv"),
      path("metrics_per_cell_exclude_peak_gene_regions_only_background.csv"),
      path("metrics_per_cell_exclude_nothing.csv"),
      path("metrics_per_cell_exclude_peak_regions_only_background.csv"),
      path("metrics_per_cell_exclude_peak_gene_regions.csv"),
      path("metrics_per_cell_only_peak_regions.csv"),
      path(rt_ranks)
    
    output:
    path "*.{pdf,png}"
    
    script:
    """
    generate_plots.py \\
        ${rt_ranks} \\
        ${metrics} \\
        ./
    """
}



// Error handling
workflow.onError {
    println "Pipeline execution stopped with error: ${workflow.errorMessage}"
}