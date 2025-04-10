# --- Data Preparation Rules ---
rule download_and_index_fasta:
    output:
        fasta = GENOME_FASTA,
        fai   = f"{GENOME_FASTA}.fai"
    params:
        url    = config["genome_fasta_url_template"].format(genome_build=config["genome_build"]),
        gz_tmp = f"data/{GENOME_BUILD}.fa.gz"
    container:
        config["utils_container"]
    resources:
        mem_mb  = RESOURCES["download_and_index_fasta"]["mem_mb"],
        runtime = RESOURCES["download_and_index_fasta"]["runtime"]
    log:
        f"{OUTPUT_DIR}/logs/data_prep/download_fasta_{config['genome_build']}.log"
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
    shell:
        """
        wget -O {output.chrom_sizes} {params.url}
        """

rule get_blacklist:
    output:
        blacklist_gz=f"data/{GENOME_BUILD}.blacklist.bed.gz",
        blacklist_bed=f"data/{GENOME_BUILD}.blacklist.bed"
    params:
        url=config["blacklist_url"].format(genome_build=GENOME_BUILD)
    shell:
        """
        wget -O {output.blacklist_gz} {params.url}
        gunzip -c {output.blacklist_gz} > {output.blacklist_bed}
        """