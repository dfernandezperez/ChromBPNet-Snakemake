# --- QC of the chrombpnet models (same steps and file names as `chrombpnet pipeline`) ---
# The full report (evaluation/<sample>_overall_report.html) is rebuilt with chrombpnet's own report
# generator once all its parts exist. Its TF-MoDISco motifs come from a random subset of the fold's
# contribution scores of all the peaks (computed anyway for the downstream analyses), instead of
# interpreting 30K peaks again on the GPU.

# Bias motifs whose marginal footprints show whether the bias was removed (as in `chrombpnet pipeline`)
BIAS_MOTIFS = os.path.join(workflow.basedir, "motifs", f"bias_motifs_{config['chromnpnet_model']['assay_type']}.tsv")

# Test set metrics (counts correlation, profile JSD) and marginal footprints of the bias motifs.
# Runs right after training, so these early QC results don't wait for the contribution scores.
rule chrombpnet_qc_gpu:
    input:
        genome          = GENOME_FASTA,
        split           = FOLD_JSON,
        model           = rules.train_chrombpnet.output.model,
        model_with_bias = rules.train_chrombpnet.output.model_with_bias,
        bigwig          = rules.train_chrombpnet.output.bigwig,
        peaks           = rules.train_chrombpnet.output.peaks,
        nonpeaks        = "preprocessing/nonpeaks/{sample}_fold_{fold}_negatives.bed",
        bias_motifs     = BIAS_MOTIFS
    output:
        metrics       = f"{MODEL_DIR}/evaluation/{{sample}}_chrombpnet_metrics.json",
        bias_response = f"{MODEL_DIR}/evaluation/{{sample}}_chrombpnet_nobias_max_bias_response.txt",
        footprints    = f"{MODEL_DIR}/auxiliary/{{sample}}_chrombpnet_nobias_footprints.h5"
    params:
        prefix_pred       = f"{MODEL_DIR}/evaluation/{{sample}}_chrombpnet",
        prefix_footprints = f"{MODEL_DIR}/evaluation/{{sample}}_chrombpnet_nobias",
        thread_env        = lambda w, threads: thread_env(threads)
    resources:
        mem_mb    = RESOURCES["chrombpnet_qc_gpu"]["mem_mb"],
        cpu       = RESOURCES["chrombpnet_qc_gpu"]["cpu"],
        runtime   = RESOURCES["chrombpnet_qc_gpu"]["runtime"],
        gres      = RESOURCES["chrombpnet_qc_gpu"]["gres"],
        slurm_partition = gpu_partition("chrombpnet_qc_gpu")
    threads:
        RESOURCES["chrombpnet_qc_gpu"]["cpu"]
    retries:
        RESOURCES["chrombpnet_qc_gpu"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/chrombpnet_qc/{sample}_fold_{fold}_gpu.log"
    shell:
        """
        export {params.thread_env}
        python -m chrombpnet.training.predict \
            -g {input.genome} \
            -b {input.bigwig} \
            -p {input.peaks} \
            -fl {input.split} \
            -m {input.model_with_bias} \
            -o {params.prefix_pred} \
            > {log} 2>&1
        chrombpnet footprints \
            -m {input.model} \
            -r {input.nonpeaks} \
            -g {input.genome} \
            -fl {input.split} \
            -pwm_f {input.bias_motifs} \
            -op {params.prefix_footprints} \
            >> {log} 2>&1
        mv {params.prefix_footprints}_footprints.h5 {output.footprints}
        """

# TF-MoDISco on a random subset of the fold's peaks (the contribution scores are stored in a random
# order, see split_contribs_regions, so the first N regions are a random subset)
rule chrombpnet_qc_modisco:
    input:
        scores  = "contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.h5",
        regions = "contribs_bw/{sample}/{sample}_fold_{fold}.interpreted_regions.bed"
    output:
        scores  = temp(f"{MODEL_DIR}/auxiliary/interpret_subsample/{{sample}}_chrombpnet_nobias.profile_scores.h5"),
        regions = f"{MODEL_DIR}/auxiliary/{{sample}}_30K_subsample_peaks.bed",
        modisco = f"{MODEL_DIR}/auxiliary/interpret_subsample/{{sample}}_modisco_results_profile_scores.h5",
        html    = f"{MODEL_DIR}/evaluation/modisco_profile/motifs.html",
        pdf     = f"{MODEL_DIR}/evaluation/{{sample}}_chrombpnet_nobias_profile.pdf"
    params:
        scripts    = SCRIPTS_DIR,
        report_dir = f"{MODEL_DIR}/evaluation/modisco_profile",
        n_regions  = config["qc_modisco"]["n_regions"],
        n_seqlets  = config["qc_modisco"]["seqlets"],
        window     = config["qc_modisco"]["window"],
        meme_db    = config["modisco_tf"]["meme_db"],
        thread_env = lambda w, threads: thread_env(threads)
    resources:
        mem_mb  = RESOURCES["chrombpnet_qc_modisco"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_qc_modisco"]["runtime"]
    threads:
        RESOURCES["chrombpnet_qc_modisco"]["cpu"]
    retries:
        RESOURCES["chrombpnet_qc_modisco"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/chrombpnet_qc/{sample}_fold_{fold}_modisco.log"
    shell:
        """
        export {params.thread_env}
        python {params.scripts}/average_contribs_h5.py --h5 {input.scores} --beds {input.regions} \
            --out-h5 {output.scores} --out-bed {output.regions} --max-rows {params.n_regions} > {log} 2>&1
        rm -rf {params.report_dir}
        modisco motifs -i {output.scores} -n {params.n_seqlets} -w {params.window} -o {output.modisco} >> {log} 2>&1
        # Default image path prefix ("./"): the overall report expects it
        modisco report -i {output.modisco} -o {params.report_dir}/ -m {params.meme_db} >> {log} 2>&1
        python -c "import chrombpnet.evaluation.modisco.convert_html_to_pdf as c; c.main('{output.html}', '{output.pdf}')" >> {log} 2>&1
        """

rule chrombpnet_report:
    input:
        flag          = rules.train_chrombpnet.output.flag,
        metrics       = rules.chrombpnet_qc_gpu.output.metrics,
        bias_response = rules.chrombpnet_qc_gpu.output.bias_response,
        modisco       = rules.chrombpnet_qc_modisco.output.html,
        pdf           = rules.chrombpnet_qc_modisco.output.pdf
    output:
        html = f"{MODEL_DIR}/evaluation/{{sample}}_overall_report.html",
        pdf  = f"{MODEL_DIR}/evaluation/{{sample}}_overall_report.pdf"
    params:
        dir   = MODEL_DIR,
        assay = config["chromnpnet_model"]["assay_type"]
    resources:
        mem_mb  = RESOURCES["chrombpnet_report"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_report"]["runtime"]
    retries:
        RESOURCES["chrombpnet_report"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/chrombpnet_qc/{sample}_fold_{fold}_report.log"
    shell:
        "python -m chrombpnet.helpers.generate_reports.make_html -id {params.dir} -d {params.assay} "
        "-fp {wildcards.sample} -icmd pipeline > {log} 2>&1"
