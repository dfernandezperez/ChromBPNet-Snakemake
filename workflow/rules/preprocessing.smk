# --- Preprocessing Rules (Per Sample) ---
rule macs3_peak_calling_pooled:
    input:
        bam = f"{OUTPUT_DIR}/preprocessing/bams/{PREFIX}_pooled.bam",
    output:
        peaks   = f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks.narrowPeak",
        summits = f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_summits.bed"
    params:
        outdir = f"{OUTPUT_DIR}/preprocessing/peaks",
        name   = f"{PREFIX}_pooled",
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
        f"{OUTPUT_DIR}/logs/macs3/{PREFIX}_pooled.log"
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

rule filter_peaks_blacklist:
    input:
        peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks.narrowPeak",
        blacklist   = rules.get_blacklist.output.blacklist_bed,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
    output:
        filtered_peaks = f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks_no_blacklist.bed"
    params:
        slop_bp = config["chromnpnet_bias"]["peak_blacklist_extension"],
        slop_tmp = f"{OUTPUT_DIR}/preprocessing/peaks/slop_tmp.bed" # Use a temporary file for the slopped blacklist
    resources:
        mem_mb  = RESOURCES["filter_peaks_blacklist"]["mem_mb"],
        runtime = RESOURCES["filter_peaks_blacklist"]["runtime"]
    container:
        config["utils_container"] # Needs bedtools
    log:
        f"{OUTPUT_DIR}/logs/filter_peaks/{PREFIX}_pooled.log"
    shell:
        """
        bedtools slop -i {input.blacklist} -g {input.chrom_sizes} -b {params.slop_bp} > {params.slop_tmp} 2>> {log}
        bedtools intersect -v -a {input.peaks} -b {params.slop_tmp} >> {log} 2>&1 > {output.filtered_peaks}
        rm {params.slop_tmp}
        """

rule chrombpnet_prep_nonpeaks:
    input:
        genome      = GENOME_FASTA,
        peaks       = rules.filter_peaks_blacklist.output.filtered_peaks, # Use filtered peaks
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        fold_json   = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        blacklist   = rules.get_blacklist.output.blacklist_bed
    output:
        nonpeaks = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{PREFIX}_fold_{{fold}}_negatives.bed",
    params:
        output_prefix = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{PREFIX}_fold_{{fold}}"
    resources:
        mem_mb  = RESOURCES["chrombpnet_prep_nonpeaks"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_prep_nonpeaks"]["runtime"]
    container:
        config["chrombpnet_container"]
    shadow:
        "minimal"
    log:
        f"{OUTPUT_DIR}/logs/chrombpnet_nonpeaks/{PREFIX}_fold_{{fold}}.log"
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