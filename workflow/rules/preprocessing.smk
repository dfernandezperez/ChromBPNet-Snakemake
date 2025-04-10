# --- Preprocessing Rules (Per Sample) ---
# Combine MACS3 peak calling for all samples first
rule macs3_peak_calling_pooled:
    input:
        bams = lambda wildcards: [config["samples"][sample] for sample in SAMPLES],
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
    container:
        config["macs3_container"]
    log:
        f"{OUTPUT_DIR}/logs/macs3/{PREFIX}_pooled.log"
    shell:
        """
        macs3 callpeak \
            -t {input.bams} \
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
        slop_bp = config["peak_blacklist_extension"]
    container:
        config["utils_container"] # Needs bedtools
    log:
        f"{OUTPUT_DIR}/logs/filter_peaks/{PREFIX}_pooled.log"
    shell:
        # Use a temporary file for the slopped blacklist
        """
        TEMP_SLOP=$(mktemp --suffix=.bed)
        bedtools slop -i {input.blacklist} -g {input.chrom_sizes} -b {params.slop_bp} > $TEMP_SLOP 2> {log}
        bedtools intersect -v -a {input.peaks} -b $TEMP_SLOP >> {log} 2>&1 > {output.filtered_peaks}
        rm $TEMP_SLOP
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
    container:
        config["chrombpnet_container"]
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