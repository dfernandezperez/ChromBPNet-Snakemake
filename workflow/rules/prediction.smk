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
         prefix   = "{sample}_fold_{fold}",
         thread_env = lambda w, threads: thread_env(threads)
     resources:
        mem_mb    = RESOURCES["pred_bw"]["mem_mb"],
        cpu       = RESOURCES["pred_bw"]["cpu"],
        runtime   = RESOURCES["pred_bw"]["runtime"],
        gres      = RESOURCES["pred_bw"]["gres"],
        slurm_partition = gpu_partition("pred_bw")
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
        export {params.thread_env}
        mkdir -p {params.out_fold}
        chrombpnet pred_bw \
        -cmb {input.model} \
        -r {input.peaks} \
        -g {input.genome} \
        -c {input.chrom_sizes} \
        -op {params.out_fold}/{params.prefix} > {log} 2>&1
        """

# --- Contribution scores (DeepSHAP) of all the peaks ---
# The slowest step (~165 regions/min per GPU): the peaks are split into `contribs.n_chunks` GPU jobs
# that can run on several GPUs at once, each running `contribs.procs_per_gpu` processes on its GPU
# (scripts/parallel_contribs.py). The chunks are then merged into the files of `chrombpnet contribs_bw`.
# The peaks are shuffled first (fixed seed, same order for all folds, so the folds can be averaged
# row by row): TF-MoDISco takes its seqlets in region order, and the chromosome-sorted peaks would
# bias them towards the first chromosomes when the seqlet limit is reached.
CONTRIB_CHUNKS = [f"{i:03d}" for i in range(config["contribs"]["n_chunks"])]
CONTRIB_CHUNK_PREFIX = "contribs_bw/{sample}/chunks/{sample}_fold_{fold}_chunk_{chunk}"

def contribs_chunk_file(path):
    return path if config["contribs"]["keep_chunks"] else temp(path)

rule split_contribs_regions:
    input:
        peaks = f"preprocessing/peaks/{{sample}}_peaks_no_blacklist_top{TOP_N_PEAKS}.bed"
    output:
        expand("contribs_bw/{{sample}}/chunks/regions_{chunk}.bed", chunk = CONTRIB_CHUNKS)
    params:
        seed = config["contribs"]["shuffle_seed"]
    localrule: True
    run:
        import numpy as np
        with open(input.peaks) as f:
            lines = [line for line in f if line.strip()]
        order = np.random.RandomState(params.seed).permutation(len(lines))
        for chunk, path in zip(np.array_split(order, len(output)), output):
            with open(path, "w") as f:
                f.writelines(lines[i] for i in chunk)

rule contribs_chunk:
    input:
        genome  = GENOME_FASTA,
        model   = "chrombpnet_models/{sample}/fold_{fold}/models/{sample}_chrombpnet_nobias.h5",
        regions = "contribs_bw/{sample}/chunks/regions_{chunk}.bed"
    output:
        regions = contribs_chunk_file(f"{CONTRIB_CHUNK_PREFIX}.interpreted_regions.bed"),
        profile = contribs_chunk_file(f"{CONTRIB_CHUNK_PREFIX}.profile_scores.h5"),
        counts  = contribs_chunk_file(f"{CONTRIB_CHUNK_PREFIX}.counts_scores.h5")
    params:
        scripts = SCRIPTS_DIR,
        prefix  = CONTRIB_CHUNK_PREFIX,
        n_procs = config["contribs"]["procs_per_gpu"]
    resources:
        mem_mb    = RESOURCES["contribs_chunk"]["mem_mb"],
        cpu       = RESOURCES["contribs_chunk"]["cpu"],
        runtime   = RESOURCES["contribs_chunk"]["runtime"],
        gres      = RESOURCES["contribs_chunk"]["gres"],
        slurm_partition = gpu_partition("contribs_chunk")
    threads:
        RESOURCES["contribs_chunk"]["cpu"]
    retries:
        RESOURCES["contribs_chunk"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/contribs_bw/{sample}_fold_{fold}_chunk_{chunk}.log"
    shell:
        """
        python {params.scripts}/parallel_contribs.py \
            --model {input.model} \
            --regions {input.regions} \
            --genome {input.genome} \
            --out-prefix {params.prefix} \
            --heads counts profile \
            --n-procs {params.n_procs} \
            --threads {threads} \
            > {log} 2>&1
        """

rule merge_contribs:
    input:
        regions     = expand(f"{CONTRIB_CHUNK_PREFIX}.interpreted_regions.bed", chunk = CONTRIB_CHUNKS, allow_missing = True),
        profile     = expand(f"{CONTRIB_CHUNK_PREFIX}.profile_scores.h5", chunk = CONTRIB_CHUNKS, allow_missing = True),
        counts      = expand(f"{CONTRIB_CHUNK_PREFIX}.counts_scores.h5", chunk = CONTRIB_CHUNKS, allow_missing = True),
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
    output:
        out_bed     = "contribs_bw/{sample}/{sample}_fold_{fold}.interpreted_regions.bed",
        profiles    = "contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.h5",
        counts      = "contribs_bw/{sample}/{sample}_fold_{fold}.counts_scores.h5",
        profiles_bw = "contribs_bw/{sample}/{sample}_fold_{fold}.profile_scores.bw",
        counts_bw   = "contribs_bw/{sample}/{sample}_fold_{fold}.counts_scores.bw"
    params:
        scripts = SCRIPTS_DIR,
        parts   = lambda w: expand(CONTRIB_CHUNK_PREFIX, chunk = CONTRIB_CHUNKS, sample = w.sample, fold = w.fold),
        prefix  = "contribs_bw/{sample}/{sample}_fold_{fold}"
    resources:
        mem_mb  = RESOURCES["merge_contribs"]["mem_mb"],
        runtime = RESOURCES["merge_contribs"]["runtime"]
    threads:
        RESOURCES["merge_contribs"]["cpu"]
    retries:
        RESOURCES["merge_contribs"]["retries"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/contribs_bw/{sample}_fold_{fold}_merge.log"
    shell:
        """
        python {params.scripts}/merge_contribs.py \
            --parts {params.parts} \
            --out-prefix {params.prefix} \
            --chrom-sizes {input.chrom_sizes} \
            --heads counts profile \
            > {log} 2>&1
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

# TF-MoDISco motifs (trimmed position frequency matrices) in MEME format, e.g. for FIMO or TOBIAS
rule modisco_meme:
     input:
         h5 = rules.modisco_tf.output.h5
     output:
         meme = "modisco_tf/{sample}/{sample}_modisco_motifs.meme"
     params:
         scripts = SCRIPTS_DIR
     resources:
        mem_mb  = RESOURCES["modisco_meme"]["mem_mb"],
        runtime = RESOURCES["modisco_meme"]["runtime"]
     retries:
        RESOURCES["modisco_meme"]["retries"]
     container:
         config["chrombpnet_container"]
     log:
         "logs/modisco_tf/{sample}_meme.log"
     shell:
        "python {params.scripts}/modisco_to_meme.py -i {input.h5} -o {output.meme} --name-prefix {wildcards.sample} > {log} 2>&1"

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
         motifs   = config["footprints"]["motifs_to_pwm"],
         thread_env = lambda w, threads: thread_env(threads)
     resources:
        mem_mb    = RESOURCES["chrombpnet_footprints"]["mem_mb"],
        cpu       = RESOURCES["chrombpnet_footprints"]["cpu"],
        runtime   = RESOURCES["chrombpnet_footprints"]["runtime"],
        gres      = RESOURCES["chrombpnet_footprints"]["gres"],
        slurm_partition = gpu_partition("chrombpnet_footprints")
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
        export {params.thread_env}
        mkdir -p {params.out_fold}
        chrombpnet footprints \
        -m {input.model} \
        -r {input.nonpeaks} \
        -g {input.genome} \
        -fl {input.split} \
        -pwm_f {params.motifs} \
        -op {params.out_fold}/{params.prefix} > {log} 2>&1
        """
