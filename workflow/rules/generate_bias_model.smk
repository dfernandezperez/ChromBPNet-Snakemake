rule merge_sort_index_bias_bams:
    input:
        bams = ALL_BIAS_BAMS
    output:
        bam = f"{OUTPUT_DIR}/bias_model/preprocessing/bams/bias_pool.bam",
        bai = f"{OUTPUT_DIR}/bias_model/preprocessing/bams/bias_pool.bam.bai"
    params:
        # Calculate the number of input bams here
        num_bams      = lambda w, input: len(input.bams),
        merged_unsort = f"{OUTPUT_DIR}/bias_model/preprocessing/bams/bias_merged_unsorted.bam"
    resources:
        mem_mb  = RESOURCES["merge_sort_index_bias_bams"]["mem_mb"],
        runtime = RESOURCES["merge_sort_index_bias_bams"]["runtime"]
    threads:
        RESOURCES["merge_sort_index_bias_bams"]["cpu"]
    container:
        config["chrombpnet_container"]
    log:
        f"{OUTPUT_DIR}/logs/bams/bias_merge_sort_index.log"
    shell:
        """
        if [ {params.num_bams} -eq 1 ]; then
            echo "Only one input BAM found, sorting and indexing..." > {log}
            samtools sort -@ {threads} -m 7G -o {output.bam} {input.bams[0]} 2>> {log}
        else
            echo "Merging {params.num_bams} BAM files..." > {log}
            # Use a temporary file for merged unsorted BAM
            # Using output.bam prefix for temp file name consistency
            MERGED_UNSORTED={params.merged_unsort}
            samtools merge -@ {threads} -f $MERGED_UNSORTED {input.bams} 2>> {log}
            echo "Sorting merged BAM..." >> {log}
            samtools sort -@ {threads} -m 7G -o {output.bam} $MERGED_UNSORTED 2>> {log}
            echo "Removing temporary merged file..." >> {log}
            rm $MERGED_UNSORTED
        fi
        echo "Indexing final BAM..." >> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        echo "BAM processing finished." >> {log}
        """

# --- Preprocessing Rules (for bias model) ---
rule macs3_bias_peak_calling:
    input:
        bam = f"{OUTPUT_DIR}/bias_model/preprocessing/bams/bias_pool.bam"
    output:
        peaks   = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_peaks.narrowPeak",
        summits = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_summits.bed"
    params:
        outdir = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks",
        name   = "bias_pool",
        gsize  = config["macs3"]["macs3_gsize"],
        pvalue = config["macs3"]["macs3_pvalue"],
        format = config["macs3"]["macs3_format"],
        model  = config["macs3"]["model_params"],
        other  = config["macs3"]["other_params"]
    threads: 
        RESOURCES["macs3_peak_calling_pooled"]["cpu"]
    resources:
        mem_mb  = RESOURCES["macs3_peak_calling_pooled"]["mem_mb"],
        runtime = RESOURCES["macs3_peak_calling_pooled"]["runtime"]
    container:
        config["macs3_container"]
    log:
        f"{OUTPUT_DIR}/logs/macs3/bias_pool.log"
    shell:
        """
        macs3 callpeak \
            -t {input.bam} \
            -f {params.format} \
            -g {params.gsize} \
            -p {params.pvalue} \
            {params.model} \
            {params.other} \
            --outdir {params.outdir} \
            -n {params.name} \
            --verbose 3 2> {log}
        """

rule filter_bias_peaks_blacklist:
    input:
        peaks       = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_peaks.narrowPeak",
        blacklist   = rules.get_blacklist.output.blacklist_bed,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
    output:
        filtered_peaks = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist.bed"
    params:
        slop_bp = config["chromnpnet_bias"]["peak_blacklist_extension"],
        slop_tmp = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_slop_tmp.bed" # Use a temporary file for the slopped blacklist
    resources:
        mem_mb  = RESOURCES["filter_peaks_blacklist"]["mem_mb"],
        runtime = RESOURCES["filter_peaks_blacklist"]["runtime"]
    container:
        config["utils_container"] # Needs bedtools
    log:
        f"{OUTPUT_DIR}/logs/filter_peaks/bias_pool.log"
    shell:
        """
        bedtools slop -i {input.blacklist} -g {input.chrom_sizes} -b {params.slop_bp} > {params.slop_tmp} 2>> {log}
        bedtools intersect -v -a {input.peaks} -b {params.slop_tmp} >> {log} 2>&1 > {output.filtered_peaks}
        rm {params.slop_tmp}
        """

rule get_top_n_peaks_bias:
    input:
        peaks = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist.bed"
    output:
        top_peaks = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist_top{TOP_N_PEAKS}.bed"
    params:
        top_n_peaks = config["macs3_top_n_peaks"],
    resources:
        mem_mb  = RESOURCES["get_top_n_peaks"]["mem_mb"],
        runtime = RESOURCES["get_top_n_peaks"]["runtime"]
    container:
        None
    log:
        f"{OUTPUT_DIR}/logs/top_n_peaks/bias_pool_top{TOP_N_PEAKS}.log"
    script:
        "../scripts/filter_top_peaks.py"

rule chrombpnet_bias_prep_nonpeaks:
    input:
        genome      = GENOME_FASTA,
        peaks       = rules.get_top_n_peaks_bias.output.top_peaks, # Use filtered and top N macs3 peaks
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        fold_json   = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        blacklist   = rules.get_blacklist.output.blacklist_bed
    output:
        nonpeaks = f"{OUTPUT_DIR}/bias_model/preprocessing/nonpeaks/bias_pool_fold_{{fold}}_negatives.bed",
    params:
        output_prefix = f"{OUTPUT_DIR}/bias_model/preprocessing/nonpeaks/bias_pool_fold_{{fold}}"
    resources:
        mem_mb  = RESOURCES["chrombpnet_prep_nonpeaks"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_prep_nonpeaks"]["runtime"]
    container:
        config["chrombpnet_container"]
    shadow:
        "minimal"
    log:
        f"{OUTPUT_DIR}/logs/chrombpnet_nonpeaks/bias_pool_fold_{{fold}}.log"
    shell:
        """
        chrombpnet prep nonpeaks \
                -g {input.genome} \
                -p {input.peaks} \
                -c {input.chrom_sizes} \
                -fl {input.fold_json} \
                -br {input.blacklist} \
                -o {params.output_prefix} \
                > {log} 2>&1
        """

# --- Model Training Rule (Per Fold) ---
rule train_bias_model:
    input:
        genome      = GENOME_FASTA,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        split       = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        bam         = f"{OUTPUT_DIR}/bias_model/preprocessing/bams/bias_pool.bam",
        nonpeaks    = f"{OUTPUT_DIR}/bias_model/preprocessing/nonpeaks/bias_pool_fold_{{fold}}_negatives.bed",
        peaks       = f"{OUTPUT_DIR}/bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
    output:
        out_dir    = directory(f"{OUTPUT_DIR}/bias_model/bias_models/fold_{{fold}}"),
        flag       = f"{OUTPUT_DIR}/bias_model/bias_models/fold_{{fold}}/complete.flag", # Add flag for rule completion
        bias_model = f"{OUTPUT_DIR}/bias_model/bias_models/fold_{{fold}}/models/bias_bias.h5",
    params:
        prefix = "pooled",
        assay  = config["chromnpnet_bias"]["assay_type"],
        extra  = config["chromnpnet_bias"]["extra_params"],
        out_dir_tmp = f"{OUTPUT_DIR}/bias_model/tmp_bias_models/fold_{{fold}}"
    resources:
        mem_mb    = RESOURCES["train_bias_model"]["mem_mb"],
        runtime   = RESOURCES["train_bias_model"]["runtime"],
        cpu       = RESOURCES["train_bias_model"]["cpu"],
        gres      = RESOURCES["train_bias_model"]["gres"],
        slurm_partition = RESOURCES["train_bias_model"]["slurm_partition"]
    retries: 
        RESOURCES["train_bias_model"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        f"{OUTPUT_DIR}/logs/train_bias/bias_fold_{{fold}}.log"
    shell:
        """
        chrombpnet bias pipeline \
                -ibam {input.bam} \
                -d {params.assay} \
                -g {input.genome} \
                -c {input.chrom_sizes} \
                -p {input.peaks} \
                -n {input.nonpeaks} \
                -fl {input.split} \
                -b 0.5 \
                -o {params.out_dir_tmp} \
                -fp {params.prefix} \
                {params.extra} \
                > {log} 2>&1
        # chrombpnet raises an error if the output folder exists at launch.
        # Snakemake creates the outdir before rule execution, so the only workaround
        # is to run the tool in a tmp folder and then move it to the actual output 
        mv {params.out_dir_tmp} {output.out_dir}
        touch {output.flag}
        """