# --- Prediction and Interpretation Rules ---

# Example: Predict on test chromosomes for a fold
rule predict_chrombpnet:
     input:
         genome=GENOME_FASTA,
         model_dir=f"{OUTPUT_DIR}/models/fold_{{fold}}/chrombpnet_model",
         splits=os.path.join(SPLITS_DIR, "fold_{fold}.json"), # Used to get test regions
         chrom_sizes=f"data/{GENOME_BUILD}.chrom.sizes"
     output:
         predictions_bw=f"{OUTPUT_DIR}/predictions/fold_{{fold}}/test_predictions.bw"
     params:
         outdir=f"{OUTPUT_DIR}/predictions/fold_{{fold}}",
         prefix=f"{PREFIX}_fold_{{fold}}_test"
     container:
         config["chrombpnet_container"]
     log:
         f"{OUTPUT_DIR}/logs/predict/fold_{{fold}}.log"
     shell:
         # Check `chrombpnet predict` syntax
         """
         chrombpnet predict -m {input.model_dir} -g {input.genome} \
                  -f {input.splits} --chrom-sizes {input.chrom_sizes} \
                  -o {params.outdir} -op {params.prefix} \
                  --predict-on test > {log} 2>&1
         # Might need to convert output HDF5 to BigWig separately if needed
         """

# Example: Interpretation
rule interpret_chrombpnet:
     input:
         genome=GENOME_FASTA,
         model_dir=f"{OUTPUT_DIR}/models/fold_{{fold}}/chrombpnet_model",
         splits=os.path.join(SPLITS_DIR, "fold_{fold}.json"), # Define regions
         peaks=f"{OUTPUT_DIR}/preprocessing/peaks/{PREFIX}_pooled_peaks.narrowPeak" # Regions to interpret
     output:
         # Interpretation often produces HDF5 or similar
         interpretation_results=f"{OUTPUT_DIR}/interpretation/fold_{{fold}}/results.h5" # Example
     params:
         outdir=f"{OUTPUT_DIR}/interpretation/fold_{{fold}}",
         prefix=f"{PREFIX}_fold_{{fold}}_interpret"
     container:
         config["chrombpnet_container"]
     log:
         f"{OUTPUT_DIR}/logs/interpret/fold_{{fold}}.log"
     shell:
          # Check `chrombpnet interpret` syntax
         """
         chrombpnet interpret -m {input.model_dir} -g {input.genome} \
                  -r {input.peaks} -f {input.splits} \
                  -o {params.outdir} -op {params.prefix} \
                  --interp-on test > {log} 2>&1
         # Add flags for specific interpretation modes (scores, modisco)
         """