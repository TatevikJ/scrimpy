#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Required parameters
params.input_paths_csv = false
params.outdir = false

// Input mode 1: RT timing data
params.rt_tsv = false

// Input mode 2: Inferring replication timing
params.rt_chrom_sizes = false
params.rt_fragments = false

// Optional parameters with defaults
params.mem_per_cpu = '40G'
params.cpus = 6

// Determine input mode
def input_mode = ""
if (params.rt_tsv) {
    input_mode = "input_rt"
} else if (params.rt_chrom_sizes && params.rt_fragments) {
    input_mode = "infer_rt"
} else {
    input_mode = "invalid"
}

// Print parameter summary
log.info """
    SCRIMPy: Single Cell Replication Inference from Multiome data using Python
    ==========================================================================
    input_paths_csv: ${params.input_paths_csv}       
    input_mode:      ${input_mode}
    ${input_mode == "input_rt" ? "rt_tsv:          ${params.rt_tsv}" : ""}
    ${input_mode == "infer_rt" ? "rt_chrom_sizes:  ${params.rt_chrom_sizes}" : ""}
    ${input_mode == "infer_rt" ? "rt_fragments:    ${params.rt_fragments}" : ""}
    outdir:          ${params.outdir}
    """

// Validate required inputs
if (!params.input_paths_csv || !params.outdir) {
    error "Missing required parameters. Please provide input_paths_csv and outdir."
}

// Validate input mode
if (input_mode == "invalid") {
    error """
    Invalid input configuration. Please provide:
        --rt_tsv <file>
    OR
        --rt_chrom_sizes <file> --rt_fragments <file>
    """
}

// Validate mutually exclusive options
if (params.rt_tsv && (params.rt_chrom_sizes || params.rt_fragments)) {
    error "Cannot specify both rt_tsv and rt_chrom_sizes/rt_fragments. Choose one input mode."
}

// Main workflow
workflow {
    sample_info_ch = Channel
                        .fromPath(params.input_paths_csv)
                        .splitCsv(header: true)
                        .map { row ->
                            // Extract the values from each row
                            def sample = row.sample
                            def fragments = row.fragments  
                            def barcodes = row.barcodes
                            
                            // Return a tuple with the values
                            tuple(sample, file(fragments), file(barcodes))
                        }    
    
    filtered_fragments = FILTER_FRAGMENTS(sample_info_ch)

    // Get rt_tsv_ch based on input mode
    if (input_mode == "input_rt") {
        rt_tsv_ch = Channel.fromPath(params.rt_tsv)
    } else if (input_mode == "infer_rt") {
        rt_chrom_sizes_ch = Channel.fromPath(params.rt_chrom_sizes)
        rt_fragments_ch = Channel.fromPath(params.rt_fragments)

        rt_tsv_ch = INFER_RT(rt_chrom_sizes_ch, rt_fragments_ch)
    }

    rt_ranks = GET_REPLICATION_RANKS(rt_tsv_ch)

    filtered_fragments_rt_ranks_combined_ch = filtered_fragments.combine(rt_ranks)

    fragment_ranks = GET_FRAGMENT_RANKS(filtered_fragments_rt_ranks_combined_ch)

    processed_fragments = PROCESS_FRAGMENT_RANKS(fragment_ranks)

    overlapping_fragments = GET_OVERLAPPING_FRAGMENTS(filtered_fragments)

    // Join by sample name
    overlapping_fragments_processed_fragments_joined_ch = overlapping_fragments.join(processed_fragments, by: 0)

    //overlapping_fragments_processed_fragments_joined_ch.view()

    metrics = PROCESS_OVERLAPPING_FRAGMENTS(overlapping_fragments_processed_fragments_joined_ch)

    metrics_combined_ch = metrics.combine(rt_ranks)
}


// Process 1: Filter atac fragments by provided list of barcodes and keep only fragments from autosomal chromosomes for each sample
process FILTER_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), 
          path(fragments), 
          path(barcodes)
    
    output:
    tuple val(sample), 
          path("atac_fragments_filtered.tsv.gz"), 
          path(barcodes)
    
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

// Process 2*: Infer RT from input fragment coverage per 50kb bins of genome
process INFER_RT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path rt_chrom_sizes
    path rt_fragments
    
    output:
    path "inferred_rt.tsv"
    
    script:
    """
    # Create genomic windows
    bedtools makewindows \\
        -g ${rt_chrom_sizes} \\
        -w 50000 \\
        > genomic_windows_50kb.bed
    
    # Calculate coverage of fragments across genomic windows
    bedtools coverage \\
        -a genomic_windows_50kb.bed \\
        -b ${rt_fragments} \\
        > fragment_coverage_per_bin.bed
    
    # Process coverage data with Python
    python3 << 'EOF'
    import pandas as pd
    
    # Read coverage data
    coverage_data = pd.read_csv("fragment_coverage_per_bin.bed", sep="\\t", header=None)
    coverage_data.columns = ['bin_chrom', 'bin_chromStart', 'bin_chromEnd', 'num_fragments', 'num_bases_covered', 'total_bases', 'coverage_fraction']
    
    coverage_data[['bin_chrom', 'bin_chromStart', 'num_fragments']].to_csv("inferred_rt.tsv",
                                                                          sep='\\t',
                                                                          index=None)
    EOF
    """
}

// Process 2: Get replication ranks per bin for each sample
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

// Process 3: Get fragment ranks for each sample
process GET_FRAGMENT_RANKS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), 
          path(filtered_fragments),
          path(barcodes),
          path(rt_ranks)
    
    output:
    tuple val(sample),
          path("atac_fragments_ranks.bed.gz"),
          path(barcodes)
    
    script:
    """
    # Get ranks for all fragments
    bedtools intersect -wo \\
        -a ${filtered_fragments} \\
        -b ${rt_ranks} \\
        > atac_fragments_ranks.bed

    gzip atac_fragments_ranks.bed
    """
}

// Process 4: Process fragment ranks: 1. filter fragment ranks, 2. get rank counts per cell
process PROCESS_FRAGMENT_RANKS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample),
          path(fragments_ranks),
          path(barcodes)
    
    output:
    tuple val(sample),
          path("atac_fragments_ranks_filtered.tsv.gz"),
          path("atac_fragments_ranks_counts_per_cell.tsv")
    script:
    """
    process_fragments_ranks.py \\
        ${fragments_ranks} \\
        ${barcodes} \\
        atac_fragments_ranks_filtered.tsv \\
        atac_fragments_ranks_counts_per_cell.tsv 

    gzip atac_fragments_ranks_filtered.tsv
    """
}


// Process 5: Get overlapping fragments for each sample
process GET_OVERLAPPING_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), 
          path(filtered_fragments),
          path(barcodes)

    output:
    tuple val(sample),
          path("atac_overlapping_fragments.bed"),
          path(filtered_fragments),
          path(barcodes)
          
    script:
    """    
    #set -euo pipefail
    #### Get overlapping fragments - exclude nothing
    #bedtools intersect \\
    #    -a ${filtered_fragments} \\
    #    -b ${filtered_fragments} \\
    #    -wo -sorted | \\
    #awk '\$1 == \$6 && \$4 == \$9 && (\$2 < \$7 || (\$2 == \$7 && \$3 < \$8))' - \\
    #    > atac_overlapping_fragment_pairs.bed
    #
    #awk '{print \$1, \$2, \$3, \$4, \$5; print \$6, \$7, \$8, \$9, \$10}' OFS="\\t" atac_overlapping_fragment_pairs.bed | \\
    #sort -k1,1 -k2,2n -k3,3n -u \\
    #    > atac_overlapping_fragments.bed



    set -euo pipefail
    mkdir -p "sort_tmp"
    trap "rm -rf sort_tmp" EXIT
    
    bedtools intersect -wo -sorted -a ${filtered_fragments} -b ${filtered_fragments} | \\
    awk '
      \$1 == \$6 && \$4 == \$9 && (\$2 < \$7 || (\$2 == \$7 && \$3 < \$8)) {
        print \$1, \$2, \$3, \$4, \$5; print \$6, \$7, \$8, \$9, \$10
      }
    ' OFS="\\t" | \\
    sort -k1,1 -k2,2n -k3,3n -u -S 8G -T sort_tmp \\
    > atac_overlapping_fragments.bed

    """
}

// Process 7: Process overlapping fragments and generate metrics
process PROCESS_OVERLAPPING_FRAGMENTS {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    
    input:
    tuple val(sample),
          path(overlapping_fragments),
          path(filtered_fragments),
          path(barcodes),
          path(fragments_ranks_filtered),
          path(ranks_counts_per_cell)

    output:
    tuple val(sample), 
          path("metrics_per_cell.csv")

    script:
    """    
    process_overlapping_fragments.py \\
        ${overlapping_fragments} \\
        ${barcodes} \\
        ${filtered_fragments} \\
        ${fragments_ranks_filtered} \\
        metrics_per_cell.csv
    """
}

// Error handling
workflow.onError {
    println "Pipeline execution stopped with error: ${workflow.errorMessage}"
}