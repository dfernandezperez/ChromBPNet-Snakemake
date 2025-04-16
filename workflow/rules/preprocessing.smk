# --- Preprocessing Rules (Per Sample) ---
rule macs3_peak_calling:
    input:
        bam = f"{OUTPUT_DIR}/preprocessing/bams/{{sample}}.bam"
    output:
        peaks   = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks.narrowPeak",
        summits = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_summits.bed"
    params:
        outdir = f"{OUTPUT_DIR}/preprocessing/peaks",
        name   = lambda w: w.sample,
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
        f"{OUTPUT_DIR}/logs/macs3/{{sample}}_pooled.log"
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
        peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks.narrowPeak",
        blacklist   = rules.get_blacklist.output.blacklist_bed,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
    output:
        filtered_peaks = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks_no_blacklist.bed"
    params:
        slop_bp = config["chromnpnet_bias"]["peak_blacklist_extension"],
        slop_tmp = lambda w: f"{OUTPUT_DIR}/preprocessing/peaks/{{w.sample}}_slop_tmp.bed" # Use a temporary file for the slopped blacklist
    resources:
        mem_mb  = RESOURCES["filter_peaks_blacklist"]["mem_mb"],
        runtime = RESOURCES["filter_peaks_blacklist"]["runtime"]
    container:
        config["utils_container"] # Needs bedtools
    log:
        f"{OUTPUT_DIR}/logs/filter_peaks/{{sample}}.log"
    shell:
        """
        bedtools slop -i {input.blacklist} -g {input.chrom_sizes} -b {params.slop_bp} > {params.slop_tmp} 2>> {log}
        bedtools intersect -v -a {input.peaks} -b {params.slop_tmp} >> {log} 2>&1 > {output.filtered_peaks}
        rm {params.slop_tmp}
        """

rule get_top_n_peaks:
    input:
        peaks = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks_no_blacklist.bed"
    output:
        top_peaks = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed" # Embed N in filename if desired
    params:
        top_n_peaks = config["macs3_top_n_peaks"],
    resources:
        mem_mb  = RESOURCES["get_top_n_peaks"]["mem_mb"],
        runtime = RESOURCES["get_top_n_peaks"]["runtime"]
    container:
        None
    log:
        f"{OUTPUT_DIR}/logs/top_n_peaks/{{sample}}_top{TOP_N_PEAKS}.log" # Match output filename
    shell:
        """
        score_threshold=$(cut -f7 {input.peaks} | sort -nr | awk -v n={params.top_n_peaks} "NR==n{print; exit} ENif(NR<n) print $1}")

        echo "Score threshold for top {params.top_n_peaks} peaks: $score_threshold" >> {log}

        awk -v "BEGIN{OFS="\t"} $7 >= $score_threshold" {input.peaks} > {output.top_peaks} 2>> {log}

        echo "Filtered peaks saved to {output.top_peaks}" >> {log}
        """

rule chrombpnet_prep_nonpeaks:
    input:
        genome      = GENOME_FASTA,
        peaks       = rules.get_top_n_peaks.output.top_peaks, # Use filtered and top N macs3 peaks
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        fold_json   = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        blacklist   = rules.get_blacklist.output.blacklist_bed
    output:
        nonpeaks = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{{sample}}_fold_{{fold}}_negatives.bed",
    params:
        output_prefix = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{{sample}}_fold_{{fold}}"
    resources:
        mem_mb  = RESOURCES["chrombpnet_prep_nonpeaks"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_prep_nonpeaks"]["runtime"]
    container:
        config["chrombpnet_container"]
    shadow:
        "minimal"
    log:
        f"{OUTPUT_DIR}/logs/chrombpnet_nonpeaks/{{sample}}_fold_{{fold}}.log"
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