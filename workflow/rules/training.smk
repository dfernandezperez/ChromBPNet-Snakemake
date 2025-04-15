# --- Model Training Rules (Per Fold) ---
rule train_bias_model:
    input:
        genome      = GENOME_FASTA,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        split       = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        bam         = f"{OUTPUT_DIR}/preprocessing/bams/{PREFIX}_pooled.bam",
        peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks_no_blacklist.bed",
        nonpeaks    = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{PREFIX}_fold_{{fold}}_negatives.bed"
    output:
        out_dir    = directory(f"{OUTPUT_DIR}/bias_models/fold_{{fold}}"),
        flag       = f"{OUTPUT_DIR}/bias_models/fold_{{fold}}/complete.flag", # Add flag for rule completion
        bias_model = f"{OUTPUT_DIR}/bias_models/fold_{{fold}}/models/{PREFIX}_bias.h5",
    params:
        prefix = f"{PREFIX}",
        assay  = config["chromnpnet_bias"]["assay_type"],
        extra  = config["chromnpnet_bias"]["extra_params"],
        out_dir_tmp = f"{OUTPUT_DIR}/tmp_bias_models/fold_{{fold}}"
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
        f"{OUTPUT_DIR}/logs/train_bias/{PREFIX}_fold_{{fold}}.log"
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

rule train_chrombpnet:
    input:
        genome      = GENOME_FASTA,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        split       = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        bam         = f"{OUTPUT_DIR}/preprocessing/bams/{PREFIX}_pooled.bam",
        peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks_no_blacklist.bed",
        nonpeaks    = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{PREFIX}_fold_{{fold}}_negatives.bed",
        bias_model  = f"{OUTPUT_DIR}/bias_models/fold_{{fold}}/models/{PREFIX}_bias.h5"
    output:
        out_dir = directory(f"{OUTPUT_DIR}/chrombpnet_models/{PREFIX}/fold_{{fold}}"),
        flag    = f"{OUTPUT_DIR}/chrombpnet_models/{PREFIX}/fold_{{fold}}/complete.flag", # Add flag for rule completion
        model   = f"{OUTPUT_DIR}/chrombpnet_models/{PREFIX}/fold_{{fold}}/models/{PREFIX}_chrombpnet_nobias.h5",
    params:
        prefix = f"{PREFIX}",
        assay  = config["chromnpnet_bias"]["assay_type"],
        extra  = config["chromnpnet_bias"]["extra_params"],
        out_dir_tmp = f"{OUTPUT_DIR}/tmp_models/{PREFIX}/fold_{{fold}}"
    resources:
        mem_mb    = RESOURCES["train_chrombpnet"]["mem_mb"],
        runtime   = RESOURCES["train_chrombpnet"]["runtime"],
        cpu       = RESOURCES["train_chrombpnet"]["cpu"],
        gres      = RESOURCES["train_chrombpnet"]["gres"],
        slurm_partition = RESOURCES["train_chrombpnet"]["slurm_partition"]
    retries: 
        RESOURCES["train_chrombpnet"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        f"{OUTPUT_DIR}/logs/train_chrombpnet/fold_{{fold}}.log"
    shell:
        """
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
                > {log} 2>&1
        # chrombpnet raises an error if the output folder exists at launch.
        # Snakemake creates the outdir before rule execution, so the only workaround
        # is to run the tool in a tmp folder and then move it to the actual output 
        mv {params.out_dir_tmp} {output.out_dir}
        touch {output.flag}
        """
