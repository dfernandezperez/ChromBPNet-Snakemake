# `chrombpnet train` only trains the model (plus its training report). The QC that `chrombpnet pipeline`
# runs afterwards is done by the rules in qc.smk, its TF-MoDISco step from the contribution scores of
# all the peaks (prediction.smk) instead of interpreting 30K peaks again on the GPU.
MODEL_DIR = "chrombpnet_models/{sample}/fold_{fold}"

rule train_chrombpnet:
    input:
        genome      = GENOME_FASTA,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        split       = FOLD_JSON,
        bam         = "preprocessing/bams/{sample}.bam",
        nonpeaks    = "preprocessing/nonpeaks/{sample}_fold_{fold}_negatives.bed",
        peaks       = f"preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
        bias_model  = "bias_model/bias_models/fold_{fold}/models/bias.h5"
    output:
        flag    = f"{MODEL_DIR}/complete.flag", # Add flag for rule completion
        model   = f"{MODEL_DIR}/models/{{sample}}_chrombpnet_nobias.h5",
        model_with_bias = f"{MODEL_DIR}/models/{{sample}}_chrombpnet.h5",
        bigwig  = f"{MODEL_DIR}/auxiliary/{{sample}}_data_unstranded.bw", # Observed Tn5 insertions, used by the QC
        peaks   = f"{MODEL_DIR}/auxiliary/{{sample}}_filtered.peaks.bed",
        report  = f"{MODEL_DIR}/evaluation/{{sample}}_train_report.html", # Training-only report (the full one is made by chrombpnet_report)
    params:
        prefix = "{sample}",
        assay  = config["chromnpnet_model"]["assay_type"],
        extra  = config["chromnpnet_model"]["extra_params"],
        out_dir = MODEL_DIR,
        out_dir_tmp = "tmp_models/{sample}/fold_{fold}",
        thread_env = lambda w, threads: thread_env(threads)
    resources:
        mem_mb    = RESOURCES["train_chrombpnet"]["mem_mb"],
        cpu       = RESOURCES["train_chrombpnet"]["cpu"],
        runtime   = RESOURCES["train_chrombpnet"]["runtime"],
        gres      = RESOURCES["train_chrombpnet"]["gres"],
        slurm_partition = gpu_partition("train_chrombpnet")
    threads:
        RESOURCES["train_chrombpnet"]["cpu"]
    retries:
        RESOURCES["train_chrombpnet"]["retries"]
    priority: 50
    container:
        config["chrombpnet_container"]
    log:
        "logs/train_chrombpnet/{sample}_fold_{fold}.log"
    shell:
        """
        export {params.thread_env}
        # Delete tmpand output folder (in case exists from a failed run)
        # to avoid chrombpnet to complain
        rm -rf {params.out_dir_tmp} > {log} 2>&1
        rm -rf {params.out_dir} >> {log} 2>&1
        chrombpnet train \
                -ibam {input.bam} \
                -d {params.assay} \
                -g {input.genome} \
                -c {input.chrom_sizes} \
                -p {input.peaks} \
                -n {input.nonpeaks} \
                -fl {input.split} \
                -b {input.bias_model} \
                -o {params.out_dir_tmp} \
                -fp {params.prefix} \
                {params.extra} \
                >> {log} 2>&1
        # The full QC report (chrombpnet_report) gets the overall_report name
        mv {params.out_dir_tmp}/evaluation/{params.prefix}_overall_report.html {params.out_dir_tmp}/evaluation/{params.prefix}_train_report.html
        mv {params.out_dir_tmp}/evaluation/{params.prefix}_overall_report.pdf {params.out_dir_tmp}/evaluation/{params.prefix}_train_report.pdf
        # chrombpnet raises an error if the output folder exists at launch.
        # Snakemake creates the outdir before rule execution, so the only workaround
        # is to run the tool in a tmp folder and then move it to the actual output
        mv {params.out_dir_tmp} {params.out_dir} >> {log} 2>&1
        touch {output.flag} >> {log} 2>&1
        """


# --- Test Rule ---
# Request logs/gpu_test/tensorflow_gpu_check_successful.flag (or ..._successful_<n>.flag to run several
# copies, e.g. to check the partition routing of the GPU jobs in logs/gpu_routing.log)
rule gpu_test_rule:
    output:
        # Create a simple flag file to mark completion
        gpu_check_flag=touch("logs/gpu_test/{name}.flag")
    wildcard_constraints:
        name = r"tensorflow_gpu_check_successful(_\d+)?"
    resources:
        mem_mb    = RESOURCES["gpu_test"]["mem_mb"],
        cpu       = RESOURCES["gpu_test"]["cpu"],
        runtime   = RESOURCES["gpu_test"]["runtime"],
        gres      = RESOURCES["gpu_test"]["gres"],
        slurm_partition = gpu_partition("gpu_test")
    threads:
        RESOURCES["gpu_test"]["cpu"]
    container:
        config["chrombpnet_container"] # Use the NVIDIA container
    log:
        # Log file will contain nvidia-smi output
        "logs/gpu_test/{name}_nvidia_smi_output.log"
    script:
        "../scripts/test_gpu_chrombpnet.py"
