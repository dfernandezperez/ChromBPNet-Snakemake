# Union peaks across all conditions to generate predictions and contributions
rule merge_condition_peaks:
    input: 
        expand(f"{OUTPUT_DIR}/contribs_bw/{{sample}}/{{sample}}_fold_{{fold}}.interpreted_regions.bed", sample = SAMPLES, fold = FOLDS)
    output: 
        f"{OUTPUT_DIR}/preprocessing/consensus_peaks/all_merged.bed"
    resources:
        mem_mb    = RESOURCES["merge_condition_peaks"]["mem_mb"],
        cpu       = RESOURCES["merge_condition_peaks"]["cpu"],
        runtime   = RESOURCES["merge_condition_peaks"]["runtime"]
    message: 
        "Merging peaks across conditions"
    container:
        config["chrombpnet_container"]
    shell:
        "cat {input} | sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct | sort -k4,4 | cut -f 1-4 > {output} "

#create header for "peaks" given via run_info
rule create_peaks_header: 
    input:
        f"{OUTPUT_DIR}/preprocessing/consensus_peaks/all_merged.bed"
    output:
        f"{OUTPUT_DIR}/tobias/peaks_header.txt"
    shell:
        """
        echo "chr\tstart\tstop\tname" > {output}
        """

#---------------------- Run TOBIAS
rule tobias_bindetect_profile:
    input: 
        motifs      = "data/motif_db.txt",
        footprints  = expand(f"{OUTPUT_DIR}/contribs_bw/{{sample}}/{{sample}}_fold_{{fold}}.profile_scores.bw", sample = SAMPLES, fold = FOLDS),
        genome      = GENOME_FASTA,
        peaks       = f"{OUTPUT_DIR}/preprocessing/consensus_peaks/all_merged.bed",
        peak_header = f"{OUTPUT_DIR}/tobias/peaks_header.txt"
    output:
        directory(f"{OUTPUT_DIR}/tobias/TFBS_profile")
    resources:
        mem_mb    = RESOURCES["tobias_bindetect"]["mem_mb"],
        cpu       = RESOURCES["tobias_bindetect"]["cpu"],
        runtime   = RESOURCES["tobias_bindetect"]["runtime"]
    threads:
        RESOURCES["tobias_bindetect"]["cpu"]
    retries: 
        RESOURCES["tobias_bindetect"]["retries"]
    conda:
         "../envs/tobias.yaml"
    log:
        f"{OUTPUT_DIR}/logs/bindetect/log"
    params:
        "--cond_names " + " ".join(SAMPLES),	#comma inserts space between elements
        config.get("bindetect", "") 
    message: 
        "Running BINDetect"
    shell:
        "export TMPDIR=/tmp;" # Force to use node tmpdir instead of the one specified in the snakefile

        "TOBIAS BINDetect --motifs {input.motifs} --signals {input.footprints} --genome {input.genome} "
        "--peaks {input.peaks} --peak_header {input.peak_header} --cores {threads} --outdir {output} {params} &>> {log}; "

        "mkdir -p " + os.path.join(OUTPUT_DIR, "overview") + ";"
        "cp " + os.path.join(OUTPUT_DIR, "TFBS", "*.txt") + " " + os.path.join(OUTPUT_DIR, "overview") + ";"	#move files to overview
        "cp " + os.path.join(OUTPUT_DIR, "TFBS", "*.xlsx") + " " + os.path.join(OUTPUT_DIR, "overview") + ";"
        "cp " + os.path.join(OUTPUT_DIR, "TFBS", "*.pdf") + " " + os.path.join(OUTPUT_DIR, "overview") + ";"

rule tobias_bindetect_count:
    input: 
        motifs      = "data/motif_db.txt",
        footprints  = expand(f"{OUTPUT_DIR}/contribs_bw/{{sample}}/{{sample}}_fold_{{fold}}.counts_scores.bw", sample = SAMPLES, fold = FOLDS),
        genome      = GENOME_FASTA,
        peaks       = f"{OUTPUT_DIR}/preprocessing/consensus_peaks/all_merged.bed",
        peak_header = f"{OUTPUT_DIR}/tobias/peaks_header.txt"
    output:
        directory(f"{OUTPUT_DIR}/tobias/TFBS_count")
    resources:
        mem_mb    = RESOURCES["tobias_bindetect"]["mem_mb"],
        cpu       = RESOURCES["tobias_bindetect"]["cpu"],
        runtime   = RESOURCES["tobias_bindetect"]["runtime"]
    threads:
        RESOURCES["tobias_bindetect"]["cpu"]
    retries: 
        RESOURCES["tobias_bindetect"]["retries"]
    conda:
         "../envs/tobias.yaml"
    log:
        f"{OUTPUT_DIR}/logs/bindetect/log"
    params:
        "--cond_names " + " ".join(SAMPLES),	#comma inserts space between elements
        config.get("bindetect", "") 
    message: 
        "Running BINDetect"
    shell:
        "export TMPDIR=/tmp;" # Force to use node tmpdir instead of the one specified in the snakefile

        "TOBIAS BINDetect --motifs {input.motifs} --signals {input.footprints} --genome {input.genome} "
        "--peaks {input.peaks} --peak_header {input.peak_header} --cores {threads} --outdir {output} {params} &>> {log}; "

        "mkdir -p " + os.path.join(OUTPUT_DIR, "overview") + ";"
        "cp " + os.path.join(OUTPUT_DIR, "TFBS", "*.txt") + " " + os.path.join(OUTPUT_DIR, "overview") + ";"	#move files to overview
        "cp " + os.path.join(OUTPUT_DIR, "TFBS", "*.xlsx") + " " + os.path.join(OUTPUT_DIR, "overview") + ";"
        "cp " + os.path.join(OUTPUT_DIR, "TFBS", "*.pdf") + " " + os.path.join(OUTPUT_DIR, "overview") + ";"
