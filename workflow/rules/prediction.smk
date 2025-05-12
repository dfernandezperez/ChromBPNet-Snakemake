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
         prefix   = f"{{sample}}_fold_{{fold}}"
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
         out_fold = f"{OUTPUT_DIR}/contribs_bw/{{sample}}",
         prefix   = f"{{sample}}_fold_{{fold}}"
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
         h5 = f"{OUTPUT_DIR}/modisco_tf/{{sample}}/{{sample}}_fold_{{fold}}_modisco.h5",
         report  = directory(f"{OUTPUT_DIR}/modisco_tf/{{sample}}/reports)"
     params:
         n_seqlets = config["modisco_tf"]["seqlets"],
         meme_db   = config["modisco_tf"]["meme_db"],
         tomtom_n_match  = config["modisco_tf"]["tomtom_n_match"]
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
        # Run modisco motif discovery
        echo "Running modisco motifs" > {log}
        modisco motifs \
        -i {input.profile_scores} \
        -n {params.n_seqlets} \
        -o {output.h5} >> {log} 2>&1
        
        # Generate modisco report
        echo "Running modisco report" >> {log}
        modisco report \
        -i {input.profile_scores} \
        -m {params.meme_db} \
        -n {params.tomtom_n_match} \
        -o {output.report} \
        -s {output.report} >> {log} 2>&1
        """