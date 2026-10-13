# --- Data Preparation Rules ---
rule download_and_index_fasta:
    output:
        fasta = GENOME_FASTA,
        fai   = f"{GENOME_FASTA}.fai"
    params:
        url    = config["genome_fasta_url_template"].format(genome_build=config["genome_build"]),
        gz_tmp = f"data/{GENOME_BUILD}.fa.gz"
    container:
        config["chrombpnet_container"] # has samtools
    resources:
        mem_mb  = RESOURCES["download_and_index_fasta"]["mem_mb"],
        runtime = RESOURCES["download_and_index_fasta"]["runtime"]
    retries:
        RESOURCES["download_and_index_fasta"]["retries"]
    log:
        f"logs/data_prep/download_fasta_{config['genome_build']}.log"
    shell:
        """
        echo "Attempting to download FASTA from {params.url}" > {log}
        wget -O {params.gz_tmp} {params.url} 2>> {log} && \
        echo "Download complete. Unzipping..." >> {log} && \
        gunzip -c {params.gz_tmp} > {output.fasta} && \
        echo "Unzip complete. Indexing with samtools..." >> {log} && \
        samtools faidx {output.fasta} 2>> {log} && \
        echo "Indexing complete. Cleaning up..." >> {log} && \
        rm {params.gz_tmp} && \
        echo "FASTA download and index finished successfully." >> {log} || \
        (echo "FASTA download/index failed. Check log: {log}" && exit 1)
        """

rule get_chrom_sizes:
    output:
        chrom_sizes=f"data/{GENOME_BUILD}.chrom.sizes"
    params:
        url=config["chrom_sizes_url"].format(genome_build=GENOME_BUILD)
    resources:
        mem_mb  = RESOURCES["get_chrom_sizes"]["mem_mb"],
        runtime = RESOURCES["get_chrom_sizes"]["runtime"]
    retries:
        RESOURCES["get_chrom_sizes"]["retries"]
    log:
        f"logs/data_prep/get_chrom_sizes_{config['genome_build']}.log"
    shell:
        """
        wget -O {output.chrom_sizes} {params.url} 2> {log}
        """

rule get_blacklist:
    input:
        BLACKLIST
    output:
        blacklist_bed=f"data/{GENOME_BUILD}.blacklist.bed"
    resources:
        mem_mb  = RESOURCES["get_blacklist"]["mem_mb"],
        runtime = RESOURCES["get_blacklist"]["runtime"]
    retries:
        RESOURCES["get_blacklist"]["retries"]
    log:
        f"logs/data_prep/get_blacklist_{config['genome_build']}.log"
    shell:
        """
        gunzip -c {input} > {output} 2> {log}
        """

rule get_motif_db:
    output:
        motif_db="data/motif_db.txt"
    params:
        url=config["motif_file_url"]
    resources:
        mem_mb  = RESOURCES["get_motif_db"]["mem_mb"],
        runtime = RESOURCES["get_motif_db"]["runtime"]
    retries:
        RESOURCES["get_motif_db"]["retries"]
    log:
        "logs/data_prep/get_motif_db.log"
    shell:
        """
        wget -O {output.motif_db} {params.url} 2> {log}
        """

# Function to get BAM paths for a sample (handles single or multiple paths)
def get_bams_for_sample(wildcards):
    paths_str = SAMPLE_INFO_DF.loc[wildcards.sample, "BamPath"]
    return paths_str.split(';') # Split semicolon-separated paths

rule merge_sort_index_bams:
    input:
        bams = get_bams_for_sample
    output:
        bam = "preprocessing/bams/{sample}.bam",
        bai = "preprocessing/bams/{sample}.bam.bai"
    params:
        scripts = SCRIPTS_DIR
    resources:
        mem_mb  = RESOURCES["merge_sort_index_bams"]["mem_mb"],
        runtime = RESOURCES["merge_sort_index_bams"]["runtime"]
    threads:
        RESOURCES["merge_sort_index_bams"]["cpu"]
    container:
        config["chrombpnet_container"]
    log:
        "logs/bams/{sample}_merge_sort_index.log"
    shell:
        # Already coordinate-sorted inputs are linked or merged without sorting again
        "bash {params.scripts}/merge_sort_index_bams.sh {threads} {output.bam} {input.bams} > {log} 2>&1"

# --- Chromosome splits (only when no pre-made `splits_dir` is given) ---
if SPLITS_DIR == "":
    rule get_canonical_chrom_sizes:
        # `chrombpnet prep splits` puts every chromosome of the chrom sizes file not used
        # for test/valid into train, so keep only canonical chromosomes
        input:
            chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
        output:
            chrom_sizes = f"data/{GENOME_BUILD}.canonical.chrom.sizes"
        params:
            regex = config["chr_splits"]["canonical_chroms_regex"]
        run:
            pattern = re.compile(params.regex)
            with open(input.chrom_sizes) as f_in, open(output.chrom_sizes, "w") as f_out:
                for line in f_in:
                    if pattern.match(line.split("\t")[0]):
                        f_out.write(line)

    rule make_chr_splits:
        input:
            chrom_sizes = rules.get_canonical_chrom_sizes.output.chrom_sizes
        output:
            fold_json = "chr_splits/fold_{fold}.json"
        params:
            output_prefix = "chr_splits/fold_{fold}",
            test_chroms   = lambda w: " ".join(config["chr_splits"]["folds"][int(w.fold)]["test"]),
            valid_chroms  = lambda w: " ".join(config["chr_splits"]["folds"][int(w.fold)]["valid"])
        resources:
            mem_mb  = RESOURCES["make_chr_splits"]["mem_mb"],
            runtime = RESOURCES["make_chr_splits"]["runtime"]
        container:
            config["chrombpnet_container"]
        log:
            "logs/data_prep/make_chr_splits_fold_{fold}.log"
        shell:
            """
            chrombpnet prep splits \
                -c {input.chrom_sizes} \
                -tcr {params.test_chroms} \
                -vcr {params.valid_chroms} \
                -op {params.output_prefix} \
                > {log} 2>&1
            """
