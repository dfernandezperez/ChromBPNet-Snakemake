# --- Average predictions and contributions across folds (per sample) ---
# Recommended by ChromBPNet: models from all folds are combined before motif discovery.
# The averaged outputs feed TF-MoDISco and TOBIAS. With a single fold they are a copy of it.

rule average_pred_bw:
    input:
        bws = expand("pred_bw/{sample}_fold_{fold}_chrombpnet_nobias.bw", fold = FOLDS, allow_missing = True)
    output:
        bw = "pred_bw/{sample}_mean_chrombpnet_nobias.bw"
    params:
        scripts = SCRIPTS_DIR
    resources:
        mem_mb  = RESOURCES["average_pred_bw"]["mem_mb"],
        runtime = RESOURCES["average_pred_bw"]["runtime"]
    retries:
        RESOURCES["average_pred_bw"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/fold_average/{sample}_pred_bw.log"
    shell:
        """
        python {params.scripts}/average_bigwigs.py -o {output.bw} {input.bws} > {log} 2>&1
        """

rule average_contribs:
    input:
        beds       = expand("contribs_bw/{sample}/{sample}_fold_{fold}.interpreted_regions.bed", fold = FOLDS, allow_missing = True),
        profile_h5 = expand("contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.h5", fold = FOLDS, allow_missing = True),
        counts_h5  = expand("contribs_bw/{sample}/{sample}_fold_{fold}.counts_scores.h5", fold = FOLDS, allow_missing = True),
        profile_bw = expand("contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.bw", fold = FOLDS, allow_missing = True),
        counts_bw  = expand("contribs_bw/{sample}/{sample}_fold_{fold}.counts_scores.bw", fold = FOLDS, allow_missing = True)
    output:
        bed        = "contribs_bw/{sample}/{sample}_mean.interpreted_regions.bed",
        profile_h5 = "contribs_bw/{sample}/{sample}_mean.profile_scores.h5",
        counts_h5  = "contribs_bw/{sample}/{sample}_mean.counts_scores.h5",
        profile_bw = "contribs_bw/{sample}/{sample}_mean.profile_scores.bw",
        counts_bw  = "contribs_bw/{sample}/{sample}_mean.counts_scores.bw"
    params:
        scripts = SCRIPTS_DIR
    resources:
        mem_mb  = RESOURCES["average_contribs"]["mem_mb"],
        runtime = RESOURCES["average_contribs"]["runtime"]
    retries:
        RESOURCES["average_contribs"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/fold_average/{sample}_contribs.log"
    shell:
        """
        python {params.scripts}/average_contribs_h5.py --h5 {input.profile_h5} --beds {input.beds} \
            --out-h5 {output.profile_h5} --out-bed {output.bed} > {log} 2>&1
        python {params.scripts}/average_contribs_h5.py --h5 {input.counts_h5} --beds {input.beds} \
            --out-h5 {output.counts_h5} >> {log} 2>&1
        python {params.scripts}/average_bigwigs.py -o {output.profile_bw} {input.profile_bw} >> {log} 2>&1
        python {params.scripts}/average_bigwigs.py -o {output.counts_bw} {input.counts_bw} >> {log} 2>&1
        """
