rule pred_bw:
     input:
         genome      = GENOME_FASTA,
         model       = "chrombpnet_models/{sample}/fold_{fold}/models/{sample}_chrombpnet_nobias.h5",
         peaks       = f"preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
         chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
     output:
         pred_bw = "pred_bw/{sample}_fold_{fold}_chrombpnet_nobias.bw",
     params:
         out_fold = "pred_bw",
         prefix   = "{sample}_fold_{fold}"
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
         "logs/pred_bw/{sample}_fold_{fold}.log"
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
         model       = "chrombpnet_models/{sample}/fold_{fold}/models/{sample}_chrombpnet_nobias.h5",
         peaks       = f"preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed",
         chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
     output:
         out_bed     = "contribs_bw/{sample}/{sample}_fold_{fold}.interpreted_regions.bed",
         profiles    = "contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.h5",
         counts      = "contribs_bw/{sample}/{sample}_fold_{fold}.counts_scores.h5",
         profiles_bw = "contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.bw",
         counts_bw   = "contribs_bw/{sample}/{sample}_fold_{fold}.counts_scores.bw"
     params:
         out_fold = "contribs_bw/{sample}",
         prefix   = "{sample}_fold_{fold}"
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
         "logs/contribs_bw/{sample}_fold_{fold}.log"
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

# Run on the contribution scores averaged across folds (see fold_average.smk)
rule modisco_tf:
     input:
         profile_scores = "contribs_bw/{sample}/{sample}_mean.profile_scores.h5"
     output:
         h5 = "modisco_tf/{sample}/{sample}_modisco.h5",
         report = directory("modisco_tf/{sample}/{sample}_reports")
     params:
         n_seqlets = config["modisco_tf"]["seqlets"],
         window    = config["modisco_tf"]["window"],
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
         "logs/modisco_tf/{sample}.log"
     shell:
        """
        # Run modisco motif discovery
        echo "Running modisco motifs" > {log}
        modisco motifs \
        -i {input.profile_scores} \
        -n {params.n_seqlets} \
        -w {params.window} \
        -o {output.h5} >> {log} 2>&1
        
        # Generate modisco report
        echo "Running modisco report" >> {log}
        modisco report \
        -i {output.h5} \
        -m {params.meme_db} \
        -n {params.tomtom_n_match} \
        -o {output.report} \
        -s {output.report} >> {log} 2>&1
        """

# Marginal footprinting: motifs inserted into the non-peak regions of the fold's test chromosomes
rule chrombpnet_footprints:
     input:
         genome   = GENOME_FASTA,
         model    = "chrombpnet_models/{sample}/fold_{fold}/models/{sample}_chrombpnet_nobias.h5",
         nonpeaks = "preprocessing/nonpeaks/{sample}_fold_{fold}_negatives.bed",
         split    = FOLD_JSON
     output:
         footprints = "footprints/{sample}/{sample}_fold_{fold}_footprints.h5"
     params:
         out_fold = "footprints/{sample}",
         prefix   = "{sample}_fold_{fold}",
         motifs   = config["footprints"]["motifs_to_pwm"]
     resources:
        mem_mb    = RESOURCES["chrombpnet_footprints"]["mem_mb"],
        cpu       = RESOURCES["chrombpnet_footprints"]["cpu"],
        runtime   = RESOURCES["chrombpnet_footprints"]["runtime"],
        gres      = RESOURCES["chrombpnet_footprints"]["gres"],
        slurm_partition = RESOURCES["chrombpnet_footprints"]["slurm_partition"]
     threads:
        RESOURCES["chrombpnet_footprints"]["cpu"]
     retries: 
        RESOURCES["chrombpnet_footprints"]["retries"]
     container:
         config["chrombpnet_container"]
     log:
         "logs/footprints/{sample}_fold_{fold}.log"
     shell:
        """
        mkdir -p {params.out_fold}
        chrombpnet footprints \
        -m {input.model} \
        -r {input.nonpeaks} \
        -g {input.genome} \
        -fl {input.split} \
        -pwm_f {params.motifs} \
        -op {params.out_fold}/{params.prefix} > {log} 2>&1
        """
