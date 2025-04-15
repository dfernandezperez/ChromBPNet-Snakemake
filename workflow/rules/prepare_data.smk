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
    resources:
        mem_mb  = RESOURCES["chrombpnet_prep_nonpeaks"]["mem_mb"],
        runtime = RESOURCES["chrombpnet_prep_nonpeaks"]["runtime"]
    retries:
        RESOURCES["get_chrom_sizes"]["retries"]
    log:
        f"{OUTPUT_DIR}/logs/data_prep/get_chrom_sizes_{config['genome_build']}.log"
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
        f"{OUTPUT_DIR}/logs/data_prep/get_blacklist_{config['genome_build']}.log"
    shell:
        """
        gunzip -c {input} > {output} 2> {log}
        """

# Function to get BAM paths for a sample (handles single or multiple paths)
def get_bams_for_sample(wildcards):
    paths_str = SAMPLE_INFO_DF.loc[wildcards.sample, "BamPath"]
    return paths_str.split(';') # Split semicolon-separated paths

rule merge_sort_index_bams:
    input:
        bams = get_bams_for_sample
    output:
        bam = f"{OUTPUT_DIR}/preprocessing/bams/{{sample}}.bam",
        bai = f"{OUTPUT_DIR}/preprocessing/bams/{{sample}}.bam.bai"
    params:
        # Calculate the number of input bams here
        num_bams      = lambda w, input: len(input.bams),
        merged_unsort = f"{OUTPUT_DIR}/preprocessing/bams/{{sample}}_merged_unsorted.bam"
    resources:
        mem_mb  = RESOURCES["merge_sort_index_bams"]["mem_mb"],
        runtime = RESOURCES["merge_sort_index_bams"]["runtime"]
    threads:
        RESOURCES["merge_sort_index_bams"]["cpu"]
    container:
        config["chrombpnet_container"]
    log:
        f"{OUTPUT_DIR}/logs/bams/{{sample}}_merge_sort_index.log"
    shell:
        # Create the output directory first
        """
        if [ {params.num_bams} -eq 1 ]; then
            echo "Only one input BAM found, sorting and indexing..." > {log}
            samtools sort -@ {threads} -o {output.bam} {input.bams[0]} 2>> {log}
        else
            echo "Merging {params.num_bams} BAM files..." > {log}
            # Use a temporary file for merged unsorted BAM
            # Using output.bam prefix for temp file name consistency
            MERGED_UNSORTED={params.merged_unsort}
            samtools merge -@ {threads} -f $MERGED_UNSORTED {input.bams} 2>> {log}
            echo "Sorting merged BAM..." >> {log}
            samtools sort -@ {threads} -o {output.bam} $MERGED_UNSORTED 2>> {log}
            echo "Removing temporary merged file..." >> {log}
            rm $MERGED_UNSORTED
        fi
        echo "Indexing final BAM..." >> {log}
        samtools index {output.bam} 2>> {log}
        echo "BAM processing finished." >> {log}
        """