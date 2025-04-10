# --- Model Training Rules (Per Fold) ---

rule train_bias_model:
    input:
        genome=GENOME_FASTA,
        # !! Adjust path based on actual output of chrombpnet prep !!
        nonpeaks=f"{OUTPUT_DIR}/preprocessing/background_regions/{PREFIX}_pooled_negatives.bed",
        # !! Adjust path based on actual output of chrombpnet prep !!
        bigwig=f"{OUTPUT_DIR}/preprocessing/signal/{PREFIX}_pooled_raw.bw",
        splits=os.path.join(SPLITS_DIR, "fold_{fold}.json") # e.g., path/to/splits/fold_0.json
    output:
        # ChromBPNet saves model as a directory
        bias_model_dir=directory(f"{OUTPUT_DIR}/models/fold_{{fold}}/bias_model")
    params:
        prefix=f"{PREFIX}_fold_{{fold}}_bias",
        epochs=config["bias_model_epochs"],
        batch_size=config["train_batch_size"],
        # Filter fraction for background regions - might need adjustment
        background_fraction=config["background_fraction"]
    container:
        config["chrombpnet_container"]
    log:
        f"{OUTPUT_DIR}/logs/train_bias/fold_{{fold}}.log"
    shell:
        # Command based on wiki - check flags for bias-only training
        # May need to specify only background_regions bed file
        """
        chrombpnet bias pipeline \
                -g {input.genome} \
                -b {input.nonpeaks} \
                -s {input.bigwig} \
                -f {input.splits} \
                -o {output.bias_model_dir} \
                -op {params.prefix} \
                -e {params.epochs} \
                -bs {params.batch_size} \
                -bf {params.background_fraction} \
                 > {log} 2>&1
        """

rule train_chrombpnet:
    input:
        genome=GENOME_FASTA,
        # !! Adjust path based on actual output of chrombpnet prep !!
        peaks=f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks.narrowPeak", # Might need a processed peak file
        nonpeaks=f"{OUTPUT_DIR}/preprocessing/background_regions/{PREFIX}_pooled_negatives.bed",
        bigwig=f"{OUTPUT_DIR}/preprocessing/signal/{PREFIX}_pooled_raw.bw",
        splits=os.path.join(SPLITS_DIR, "fold_{fold}.json"),
        bias_model_dir=f"{OUTPUT_DIR}/models/fold_{{fold}}/bias_model" # Depends on previous rule
    output:
        # ChromBPNet saves model as a directory
        model_dir=directory(f"{OUTPUT_DIR}/models/fold_{{fold}}/chrombpnet_model")
    params:
        prefix=f"{PREFIX}_fold_{{fold}}_chrombpnet",
        epochs=config["train_epochs"],
        batch_size=config["train_batch_size"]
    container:
        config["chrombpnet_container"]
    log:
        f"{OUTPUT_DIR}/logs/train_chrombpnet/fold_{{fold}}.log"
    shell:
        """
        chrombpnet pipeline \
                -g {input.genome} \
                -p {input.peaks} \
                -n {input.nonpeaks} \
                -s {input.bigwig} \
                -f {input.splits} \
                -o {output.model_dir} \
                -op {params.prefix} \
                -e {params.epochs} \
                -bs {params.batch_size} \
                --bias-model-path {input.bias_model_dir} \
                 > {log} 2>&1
        """