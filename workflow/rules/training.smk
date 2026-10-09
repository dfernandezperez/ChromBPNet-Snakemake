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
        flag    = "chrombpnet_models/{sample}/fold_{fold}/complete.flag", # Add flag for rule completion
        model   = "chrombpnet_models/{sample}/fold_{fold}/models/{sample}_chrombpnet_nobias.h5",
        report  = "chrombpnet_models/{sample}/fold_{fold}/evaluation/{sample}_overall_report.html",
    params:
        prefix = "{sample}",
        assay  = config["chromnpnet_model"]["assay_type"],
        extra  = config["chromnpnet_model"]["extra_params"],
        out_dir = "chrombpnet_models/{sample}/fold_{fold}",
        out_dir_tmp = "tmp_models/{sample}/fold_{fold}"
    resources:
        mem_mb    = RESOURCES["train_chrombpnet"]["mem_mb"],
        cpu       = RESOURCES["train_chrombpnet"]["cpu"],
        runtime   = RESOURCES["train_chrombpnet"]["runtime"],
        gres      = RESOURCES["train_chrombpnet"]["gres"],
        slurm_partition = RESOURCES["train_chrombpnet"]["slurm_partition"]
    threads:
        RESOURCES["train_chrombpnet"]["cpu"]
    retries: 
        RESOURCES["train_chrombpnet"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/train_chrombpnet/{sample}_fold_{fold}.log"
    shell:
        """
        # Delete tmpand output folder (in case exists from a failed run) 
        # to avoid chrombpnet to complain
        rm -rf {params.out_dir_tmp} > {log} 2>&1
        rm -rf {params.out_dir} >> {log} 2>&1
        chrombpnet pipeline \
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
        # chrombpnet raises an error if the output folder exists at launch.
        # Snakemake creates the outdir before rule execution, so the only workaround
        # is to run the tool in a tmp folder and then move it to the actual output 
        mv {params.out_dir_tmp} {params.out_dir} >> {log} 2>&1
        touch {output.flag} >> {log} 2>&1
        """


# --- Test Rule ---
rule gpu_test_rule:
    output:
        # Create a simple flag file to mark completion
        gpu_check_flag=touch("logs/gpu_test/tensorflow_gpu_check_successful.flag")
    resources:
        mem_mb    = RESOURCES["train_chrombpnet"]["mem_mb"],
        cpu       = RESOURCES["train_chrombpnet"]["cpu"],
        runtime   = RESOURCES["train_chrombpnet"]["runtime"],
        gres      = RESOURCES["train_chrombpnet"]["gres"],
        slurm_partition = RESOURCES["train_chrombpnet"]["slurm_partition"]
    container:
        config["chrombpnet_container"] # Use the NVIDIA container
    log:
        # Log file will contain nvidia-smi output
        "logs/gpu_test/nvidia_smi_output.log"
    script:
        "../scripts/test_gpu_chrombpnet.py"