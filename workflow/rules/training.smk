rule train_chrombpnet:
    input:
        genome      = GENOME_FASTA,
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes,
        split       = os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        bam         = f"{OUTPUT_DIR}/preprocessing/bams/{{sample}}.bam",
        nonpeaks    = f"{OUTPUT_DIR}/preprocessing/nonpeaks/{{sample}}_fold_{{fold}}_negatives.bed",
        peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
        bias_model  = f"{OUTPUT_DIR}/bias_model/bias_models/fold_{{fold}}/models/bias_bias.h5"
    output:
        out_dir = directory(f"{OUTPUT_DIR}/chrombpnet_models/{{sample}}/fold_{{fold}}"),
        flag    = f"{OUTPUT_DIR}/chrombpnet_models/{{sample}}/fold_{{fold}}/complete.flag", # Add flag for rule completion
        model   = f"{OUTPUT_DIR}/chrombpnet_models/{{sample}}/fold_{{fold}}/models/{{sample}}_chrombpnet_nobias.h5",
    params:
        prefix = f"{{sample}}",
        assay  = config["chromnpnet_bias"]["assay_type"],
        extra  = config["chromnpnet_bias"]["extra_params"],
        out_dir_tmp = f"{OUTPUT_DIR}/tmp_models/{{sample}}/fold_{{fold}}"
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
        f"{OUTPUT_DIR}/logs/train_chrombpnet/{{sample}}_fold_{{fold}}.log"
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
