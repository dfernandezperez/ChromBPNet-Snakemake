rule merge_sort_index_bias_bams:
    input:
        bams = ALL_BIAS_BAMS
    output:
        bam = "bias_model/preprocessing/bams/bias_pool.bam",
        bai = "bias_model/preprocessing/bams/bias_pool.bam.bai"
    params:
        scripts = SCRIPTS_DIR
    resources:
        mem_mb  = RESOURCES["merge_sort_index_bias_bams"]["mem_mb"],
        runtime = RESOURCES["merge_sort_index_bias_bams"]["runtime"]
    threads:
        RESOURCES["merge_sort_index_bias_bams"]["cpu"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/bams/bias_merge_sort_index.log"
    shell:
        "bash {params.scripts}/merge_sort_index_bams.sh {threads} {output.bam} {input.bams} > {log} 2>&1"

# --- Preprocessing Rules (for bias model) ---
rule macs3_bias_peak_calling:
    input:
        bam = "bias_model/preprocessing/bams/bias_pool.bam"
    output:
        peaks   = "bias_model/preprocessing/peaks/bias_pool_peaks.narrowPeak",
        summits = "bias_model/preprocessing/peaks/bias_pool_summits.bed"
    params:
        outdir = "bias_model/preprocessing/peaks",
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
        "logs/macs3/bias_pool.log"
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
        peaks       = "bias_model/preprocessing/peaks/bias_pool_peaks.narrowPeak",
        blacklist   = rules.get_blacklist.output.blacklist_bed,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
    output:
        filtered_peaks = "bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist.bed"
    params:
        slop_bp = config["chromnpnet_bias"]["peak_blacklist_extension"],
        slop_tmp = "bias_model/preprocessing/peaks/bias_slop_tmp.bed" # Use a temporary file for the slopped blacklist
    resources:
        mem_mb  = RESOURCES["filter_peaks_blacklist"]["mem_mb"],
        runtime = RESOURCES["filter_peaks_blacklist"]["runtime"]
    container:
        config["utils_container"] # Needs bedtools
    log:
        "logs/filter_peaks/bias_pool.log"
    shell:
        """
        bedtools slop -i {input.blacklist} -g {input.chrom_sizes} -b {params.slop_bp} > {params.slop_tmp} 2>> {log}
        bedtools intersect -v -a {input.peaks} -b {params.slop_tmp} >> {log} 2>&1 > {output.filtered_peaks}
        rm {params.slop_tmp}
        """

rule get_top_n_peaks_bias:
    input:
        peaks = "bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist.bed"
    output:
        top_peaks = f"bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist_top{TOP_N_PEAKS}.bed"
    params:
        top_n_peaks = config["macs3_top_n_peaks"],
        canonical_chroms_regex = config["chr_splits"]["canonical_chroms_regex"]
    resources:
        mem_mb  = RESOURCES["get_top_n_peaks"]["mem_mb"],
        runtime = RESOURCES["get_top_n_peaks"]["runtime"]
    container:
        None
    log:
        f"logs/top_n_peaks/bias_pool_top{TOP_N_PEAKS}.log"
    script:
        "../scripts/filter_top_peaks.py"

rule chrombpnet_bias_prep_nonpeaks:
    input:
        genome      = GENOME_FASTA,
        peaks       = rules.get_top_n_peaks_bias.output.top_peaks, # Use filtered and top N macs3 peaks
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        fold_json   = FOLD_JSON,
        blacklist   = rules.get_blacklist.output.blacklist_bed
    output:
        nonpeaks = "bias_model/preprocessing/nonpeaks/bias_pool_fold_{fold}_negatives.bed",
    params:
        output_prefix = "bias_model/preprocessing/nonpeaks/bias_pool_fold_{fold}"
    resources:
        mem_mb  = RESOURCES["chrombpnet_prep_nonpeaks"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_prep_nonpeaks"]["runtime"]
    container:
        config["chrombpnet_container"]
    shadow:
        "minimal"
    log:
        "logs/chrombpnet_nonpeaks/bias_pool_fold_{fold}.log"
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

# --- Bias model training (per fold) ---
# `chrombpnet bias train` only trains the model, so the chrombpnet models can start training with it
# right away. Its QC (`chrombpnet bias pipeline` would run it in the same job) is split below into a GPU
# job (predictions + contribution scores of 30K peaks) and CPU jobs (TF-MoDISco + report).
BIAS_DIR = "bias_model/bias_models/fold_{fold}"

rule train_bias_model:
    input:
        genome      = GENOME_FASTA,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        split       = FOLD_JSON,
        bam         = "bias_model/preprocessing/bams/bias_pool.bam",
        nonpeaks    = "bias_model/preprocessing/nonpeaks/bias_pool_fold_{fold}_negatives.bed",
        peaks       = f"bias_model/preprocessing/peaks/bias_pool_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
    output:
        flag       = f"{BIAS_DIR}/complete.flag", # Add flag for rule completion
        bias_model = f"{BIAS_DIR}/models/bias.h5",
        bigwig     = f"{BIAS_DIR}/auxiliary/data_unstranded.bw", # Observed Tn5 insertions, used by the QC
        peaks      = f"{BIAS_DIR}/auxiliary/filtered.bias_peaks.bed",
        nonpeaks   = f"{BIAS_DIR}/auxiliary/filtered.bias_nonpeaks.bed",
        report     = f"{BIAS_DIR}/evaluation/train_report.html", # Training-only report (the full one is made by bias_report)
    params:
        assay  = config["chromnpnet_bias"]["assay_type"],
        extra  = config["chromnpnet_bias"]["extra_params"],
        out_dir = BIAS_DIR,
        out_dir_tmp = "bias_model/tmp_bias_models/fold_{fold}",
        bias_thresh_factor = config["chromnpnet_bias"]["bias_thresh_factor"],
        thread_env = lambda w, threads: thread_env(threads)
    resources:
        mem_mb    = RESOURCES["train_bias_model"]["mem_mb"],
        runtime   = RESOURCES["train_bias_model"]["runtime"],
        cpu       = RESOURCES["train_bias_model"]["cpu"],
        gres      = RESOURCES["train_bias_model"]["gres"],
        slurm_partition = gpu_partition("train_bias_model")
    threads:
        RESOURCES["train_bias_model"]["cpu"]
    retries: 
        RESOURCES["train_bias_model"]["retries"]
    priority: 100 # Every chrombpnet model needs it
    container:
        config["chrombpnet_container"]
    log:
        "logs/train_bias/bias_fold_{fold}.log"
    shell:
        """
        export {params.thread_env}
        # Delete tmp and output folder (in case exists from a failed run) 
        # to avoid chrombpnet to complain
        rm -rf {params.out_dir_tmp} > {log} 2>&1
        rm -rf {params.out_dir} >> {log} 2>&1
        chrombpnet bias train \
                -ibam {input.bam} \
                -d {params.assay} \
                -g {input.genome} \
                -c {input.chrom_sizes} \
                -p {input.peaks} \
                -n {input.nonpeaks} \
                -fl {input.split} \
                -b {params.bias_thresh_factor} \
                -o {params.out_dir_tmp} \
                {params.extra} \
                >> {log} 2>&1
        # The full QC report (bias_report) gets the overall_report name
        mv {params.out_dir_tmp}/evaluation/overall_report.html {params.out_dir_tmp}/evaluation/train_report.html
        mv {params.out_dir_tmp}/evaluation/overall_report.pdf {params.out_dir_tmp}/evaluation/train_report.pdf
        echo "bias model finished! Moving output files to output directory" >> {log} 2>&1
        # chrombpnet raises an error if the output folder exists at launch.
        # Snakemake creates the outdir before rule execution, so the only workaround
        # is to run the tool in a tmp folder and then move it to the actual output
        mv {params.out_dir_tmp} {params.out_dir} >> {log} 2>&1
        touch {output.flag} >> {log} 2>&1
        """

# --- Bias model QC (same steps and file names as `chrombpnet bias pipeline`) ---
# Predictions on the bias training peaks/nonpeaks (metrics) and contribution scores of 30K random peaks
rule bias_qc_gpu:
    input:
        genome     = GENOME_FASTA,
        split      = FOLD_JSON,
        bias_model = rules.train_bias_model.output.bias_model,
        bigwig     = rules.train_bias_model.output.bigwig,
        peaks      = rules.train_bias_model.output.peaks,
        nonpeaks   = rules.train_bias_model.output.nonpeaks
    output:
        metrics   = f"{BIAS_DIR}/evaluation/bias_metrics.json",
        subsample = f"{BIAS_DIR}/auxiliary/30K_subsample_peaks.bed",
        regions   = f"{BIAS_DIR}/auxiliary/interpret_subsample/bias.interpreted_regions.bed",
        profile   = f"{BIAS_DIR}/auxiliary/interpret_subsample/bias.profile_scores.h5",
        counts    = f"{BIAS_DIR}/auxiliary/interpret_subsample/bias.counts_scores.h5"
    params:
        scripts         = SCRIPTS_DIR,
        prefix_pred     = f"{BIAS_DIR}/evaluation/bias",
        prefix_contribs = f"{BIAS_DIR}/auxiliary/interpret_subsample/bias",
        n_regions       = config["qc_modisco"]["n_regions"],
        n_procs         = config["contribs"]["procs_per_gpu"],
        thread_env      = lambda w, threads: thread_env(threads)
    resources:
        mem_mb    = RESOURCES["bias_qc_gpu"]["mem_mb"],
        runtime   = RESOURCES["bias_qc_gpu"]["runtime"],
        cpu       = RESOURCES["bias_qc_gpu"]["cpu"],
        gres      = RESOURCES["bias_qc_gpu"]["gres"],
        slurm_partition = gpu_partition("bias_qc_gpu")
    threads:
        RESOURCES["bias_qc_gpu"]["cpu"]
    retries: 
        RESOURCES["bias_qc_gpu"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/bias_qc/bias_fold_{fold}_gpu.log"
    shell:
        """
        export {params.thread_env}
        python -m chrombpnet.training.predict \
            -g {input.genome} \
            -b {input.bigwig} \
            -p {input.peaks} \
            -n {input.nonpeaks} \
            -fl {input.split} \
            -m {input.bias_model} \
            -o {params.prefix_pred} \
            > {log} 2>&1
        python {params.scripts}/subsample_bed.py {input.peaks} {output.subsample} --n {params.n_regions} --seed 1234 >> {log} 2>&1
        python {params.scripts}/parallel_contribs.py \
            --model {input.bias_model} \
            --regions {output.subsample} \
            --genome {input.genome} \
            --out-prefix {params.prefix_contribs} \
            --heads counts profile \
            --n-procs {params.n_procs} \
            --threads {threads} \
            >> {log} 2>&1
        """

# TF-MoDISco on the bias model contributions: it should only find Tn5 (or DNase) bias motifs
rule bias_qc_modisco:
    input:
        scores  = f"{BIAS_DIR}/auxiliary/interpret_subsample/bias.{{head}}_scores.h5",
        meme_db = rules.get_modisco_meme_db.output.meme_db
    output:
        modisco = f"{BIAS_DIR}/auxiliary/interpret_subsample/modisco_results_{{head}}_scores.h5",
        html    = f"{BIAS_DIR}/evaluation/modisco_{{head}}/motifs.html"
    params:
        report_dir   = f"{BIAS_DIR}/evaluation/modisco_{{head}}",
        n_seqlets    = config["qc_modisco"]["seqlets"],
        window       = config["qc_modisco"]["window"],
        motifs_args  = MODISCO_MOTIFS_ARGS,
        thread_env   = lambda w, threads: thread_env(threads)
    resources:
        mem_mb  = RESOURCES["bias_qc_modisco"]["mem_mb"],
        runtime = RESOURCES["bias_qc_modisco"]["runtime"]
    threads:
        RESOURCES["bias_qc_modisco"]["cpu"]
    retries:
        RESOURCES["bias_qc_modisco"]["retries"]
    conda:
        MODISCO_ENV
    log:
        "logs/bias_qc/bias_fold_{fold}_modisco_{head}.log"
    shell:
        """
        export {params.thread_env}
        rm -rf {params.report_dir}
        modisco motifs -i {input.scores} -n {params.n_seqlets} -w {params.window} {params.motifs_args} -o {output.modisco} > {log} 2>&1
        # Default image path prefix ("./"): the overall report expects it
        modisco report -i {output.modisco} -o {params.report_dir}/ -m {input.meme_db} >> {log} 2>&1
        """

rule bias_report:
    input:
        flag     = rules.train_bias_model.output.flag,
        metrics  = rules.bias_qc_gpu.output.metrics,
        modisco  = expand(f"{BIAS_DIR}/evaluation/modisco_{{head}}/motifs.html", head = ["profile", "counts"], allow_missing = True)
    output:
        html = f"{BIAS_DIR}/evaluation/overall_report.html", # Check it to tune bias_thresh_factor
        pdf  = f"{BIAS_DIR}/evaluation/overall_report.pdf",
        # PDFs of the TF-MoDISco reports, as `chrombpnet bias pipeline` writes them
        modisco_pdfs = expand(f"{BIAS_DIR}/evaluation/bias_{{head}}.pdf", head = ["profile", "counts"], allow_missing = True)
    params:
        dir = BIAS_DIR
    resources:
        mem_mb  = RESOURCES["bias_report"]["mem_mb"],
        runtime = RESOURCES["bias_report"]["runtime"]
    retries: 
        RESOURCES["bias_report"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/bias_qc/bias_fold_{fold}_report.log"
    shell:
        """
        echo "Converting the TF-MoDISco reports to PDF" > {log}
        for head in profile counts; do
            python -c "import chrombpnet.evaluation.modisco.convert_html_to_pdf as c; c.main('{params.dir}/evaluation/modisco_$head/motifs.html', '{params.dir}/evaluation/bias_$head.pdf')" >> {log} 2>&1
        done
        python -m chrombpnet.helpers.generate_reports.make_html_bias -id {params.dir} -icmd pipeline >> {log} 2>&1
        """
