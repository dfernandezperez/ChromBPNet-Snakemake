rule pred_bw:
     input:
         genome      = GENOME_FASTA,
         model       = f"{OUTPUT_DIR}/chrombpnet_models/{{sample}}/fold_{{fold}}/models/{{sample}}_chrombpnet_nobias.h5",
         peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
         chrom_sizes = f"data/{GENOME_BUILD}.chrom.sizes"
     output:
         pred_bw = f"{OUTPUT_DIR}/pred_bw/{{sample}}_fold_{{fold}}_chrombpnet_nobias.bw",
     params:
         out_fold = f"{OUTPUT_DIR}/pred_bw",
         prefix   = lambda w: f"{{w.sample}}_fold_{{w.fold}}"
     resources:
        mem_mb    = RESOURCES["pred_bw"]["mem_mb"],
        cpu       = RESOURCES["pred_bw"]["cpu"],
        runtime   = RESOURCES["pred_bw"]["runtime"],
        gres      = RESOURCES["pred_bw"]["gres"],
        slurm_partition = RESOURCES["pred_bw"]["slurm_partition"]
     threads:
        RESOURCES["pred_bw"]["cpu"]
     retries: 
        RESOURCES["pred_bw"]["retries"]
     container:
         config["chrombpnet_container"]
     log:
         f"{OUTPUT_DIR}/logs/pred_bw/{{sample}}_fold_{{fold}}.log"
     shell:
        """
        mkdir -p {params.out_fold}
        chrombpnet pred_bw \
        -cmb {input.model} \
        -r {input.peaks} \
        -g {input.genome} \
        -c {input.chrom_sizes} \
        -op {params.out_fold}/{params.prefix} > {log} 2>&1
        """

rule contribs_bw:
     input:
         genome      = GENOME_FASTA,
         model       = f"{OUTPUT_DIR}/chrombpnet_models/{{sample}}/fold_{{fold}}/models/{{sample}}_chrombpnet_nobias.h5",
         peaks       = f"{OUTPUT_DIR}/preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
         chrom_sizes = f"data/{GENOME_BUILD}.chrom.sizes"
     output:
         profiles = f"{OUTPUT_DIR}/contribs_bw/{{sample}}/{{sample}}_fold_{{fold}}.profile_scores.h5",
         counts   = f"{OUTPUT_DIR}/contribs_bw/{{sample}}/{{sample}}_fold_{{fold}}.counts_scores.h5"
     params:
         out_fold = f"{OUTPUT_DIR}/contribs_b/{{sample}}/",
         prefix   = lambda w: f"{{w.sample}}_fold_{{w.fold}}"
     resources:
        mem_mb    = RESOURCES["contribs_bw"]["mem_mb"],
        cpu       = RESOURCES["contribs_bw"]["cpu"],
        runtime   = RESOURCES["contribs_bw"]["runtime"],
        gres      = RESOURCES["contribs_bw"]["gres"],
        slurm_partition = RESOURCES["contribs_bw"]["slurm_partition"]
     threads:
        RESOURCES["contribs_bw"]["cpu"]
     retries: 
        RESOURCES["contribs_bw"]["retries"]
     container:
         config["chrombpnet_container"]
     log:
         f"{OUTPUT_DIR}/logs/contribs_bw/{{sample}}_fold_{{fold}}.log"
     shell:
        """
        mkdir -p {params.out_fold}
        chrombpnet contribs_bw \
        -m {input.model} \
        -r {input.peaks} \
        -g {input.genome} \
        -c {input.chrom_sizes} \
        -op {params.out_fold}/{params.prefix} > {log} 2>&1
        """

rule modisco_tf:
     input:
         profile_scores = f"{OUTPUT_DIR}/contribs_bw/{{sample}}/{{sample}}_fold_{{fold}}.profile_scores.h5"
     output:
         modisco = f"{OUTPUT_DIR}/modisco_tf/{{sample}}/{{sample}}_fold_{{fold}}_modisco.h5"
     params:
         out_fold = f"{OUTPUT_DIR}/modisco_tf/{{sample}}/",
         prefix   = lambda w: f"{{w.sample}}_fold_{{w.fold}}"
     resources:
        mem_mb    = RESOURCES["modisco_tf"]["mem_mb"],
        cpu       = RESOURCES["modisco_tf"]["cpu"],
        runtime   = RESOURCES["modisco_tf"]["runtime"]
     threads:
        RESOURCES["modisco_tf"]["cpu"]
     retries: 
        RESOURCES["modisco_tf"]["retries"]
     container:
         config["chrombpnet_container"]
     log:
         f"{OUTPUT_DIR}/logs/modisco_tf/{{sample}}_fold_{{fold}}.log"
     shell:
        """
        mkdir -p {params.out_fold}
        modisco motifs \
        -i {input.profile_scores} \
        -n 100000 \
        -op {params.out_fold}/{params.prefix} > {log} 2>&1
        """