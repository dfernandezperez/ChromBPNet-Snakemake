# Union peaks across all conditions to generate predictions and contributions
rule merge_condition_peaks:
    input: 
        expand("contribs_bw/{sample}/{sample}_mean.interpreted_regions.bed", sample = SAMPLES)
    output: 
        "preprocessing/consensus_peaks/all_merged.bed"
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
        "preprocessing/consensus_peaks/all_merged.bed"
    output:
        "tobias/peaks_header.txt"
    shell:
        """
        echo "chr\tstart\tstop\tname" > {output}
        """

#---------------------- Run TOBIAS
rule tobias_bindetect_profile:
    input: 
        motifs      = "data/motif_db.txt",
        footprints  = expand("contribs_bw/{sample}/{sample}_mean.profile_scores.bw", sample = SAMPLES), # One fold-averaged signal per condition
        genome      = GENOME_FASTA,
        peaks       = "preprocessing/consensus_peaks/all_merged.bed",
        peak_header = "tobias/peaks_header.txt"
    output:
        directory("tobias/TFBS_profile")
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
        "logs/bindetect/profile.log"
    params:
        "--cond_names " + " ".join(SAMPLES), # One condition name per signal, in the same order
        config.get("bindetect", "") 
    message: 
        "Running BINDetect"
    shell:
        "export TMPDIR=/tmp;" # Force to use node tmpdir instead of the one specified in the snakefile

        "TOBIAS BINDetect --motifs {input.motifs} --signals {input.footprints} --genome {input.genome} "
        "--peaks {input.peaks} --peak_header {input.peak_header} --cores {threads} --outdir {output} {params} &>> {log}; "

        "mkdir -p tobias/overview_profile;"
        "cp tobias/TFBS_profile/*.txt tobias/overview_profile;"	#move files to overview
        "cp tobias/TFBS_profile/*.xlsx tobias/overview_profile;"
        "cp tobias/TFBS_profile/*.pdf tobias/overview_profile;"

rule tobias_bindetect_count:
    input: 
        motifs      = "data/motif_db.txt",
        footprints  = expand("contribs_bw/{sample}/{sample}_mean.counts_scores.bw", sample = SAMPLES), # One fold-averaged signal per condition
        genome      = GENOME_FASTA,
        peaks       = "preprocessing/consensus_peaks/all_merged.bed",
        peak_header = "tobias/peaks_header.txt"
    output:
        directory("tobias/TFBS_count")
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
        "logs/bindetect/counts.log"
    params:
        "--cond_names " + " ".join(SAMPLES), # One condition name per signal, in the same order
        config.get("bindetect", "") 
    message: 
        "Running BINDetect"
    shell:
        "export TMPDIR=/tmp;" # Force to use node tmpdir instead of the one specified in the snakefile

        "TOBIAS BINDetect --motifs {input.motifs} --signals {input.footprints} --genome {input.genome} "
        "--peaks {input.peaks} --peak_header {input.peak_header} --cores {threads} --outdir {output} {params} &>> {log}; "

        "mkdir -p tobias/overview_counts;"
        "cp tobias/TFBS_count/*.txt tobias/overview_counts;"	#move files to overview
        "cp tobias/TFBS_count/*.xlsx tobias/overview_counts;"
        "cp tobias/TFBS_count/*.pdf tobias/overview_counts;"