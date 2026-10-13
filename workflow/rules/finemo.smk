# --- Motif instances (Fi-NeMo, https://github.com/kundajelab/Fi-NeMo) ---
# Hits of each sample's TF-MoDISco motifs in its fold-averaged profile contribution scores (the same
# scores TF-MoDISco used). Main outputs:
#   finemo/<sample>/hits/hits.tsv (all hits), hits_unique.tsv (overlapping peaks deduplicated), hits.bed
#   finemo/<sample>/report/motif_occurrences.tsv (hits of each motif per peak), report.html
rule finemo_extract_regions:
    input:
        scores  = "contribs_bw/{sample}/{sample}_mean.profile_scores.h5",
        regions = "contribs_bw/{sample}/{sample}_mean.interpreted_regions.bed" # Same row order as the h5
    output:
        regions = "finemo/{sample}/{sample}_regions.npz"
    params:
        width = config["finemo"]["region_width"]
    resources:
        mem_mb  = RESOURCES["finemo_extract_regions"]["mem_mb"],
        runtime = RESOURCES["finemo_extract_regions"]["runtime"]
    retries:
        RESOURCES["finemo_extract_regions"]["retries"]
    conda:
        "../envs/finemo.yaml"
    log:
        "logs/finemo/{sample}_extract_regions.log"
    shell:
        "finemo extract-regions-chrombpnet-h5 -c {input.scores} -p {input.regions} -o {output.regions} "
        "-w {params.width} > {log} 2>&1"

rule finemo_call_hits:
    input:
        regions = rules.finemo_extract_regions.output.regions,
        modisco = "modisco_tf/{sample}/{sample}_modisco.h5"
    output:
        hits        = "finemo/{sample}/hits/hits.tsv",
        hits_unique = "finemo/{sample}/hits/hits_unique.tsv",
        hits_bed    = "finemo/{sample}/hits/hits.bed",
        motifs      = "finemo/{sample}/hits/motif_data.tsv"
    params:
        out_dir = "finemo/{sample}/hits",
        extra   = config["finemo"]["call_hits_params"]
    resources:
        mem_mb    = RESOURCES["finemo_call_hits"]["mem_mb"],
        cpu       = RESOURCES["finemo_call_hits"]["cpu"],
        runtime   = RESOURCES["finemo_call_hits"]["runtime"],
        gres      = RESOURCES["finemo_call_hits"]["gres"],
        slurm_partition = gpu_partition("finemo_call_hits")
    threads:
        RESOURCES["finemo_call_hits"]["cpu"]
    retries:
        RESOURCES["finemo_call_hits"]["retries"]
    conda:
        "../envs/finemo.yaml"
    log:
        "logs/finemo/{sample}_call_hits.log"
    shell:
        "finemo call-hits -r {input.regions} -m {input.modisco} -o {params.out_dir} {params.extra} > {log} 2>&1"

rule finemo_report:
    input:
        regions = rules.finemo_extract_regions.output.regions,
        hits    = rules.finemo_call_hits.output.hits,
        modisco = "modisco_tf/{sample}/{sample}_modisco.h5"
    output:
        report      = "finemo/{sample}/report/report.html",
        occurrences = "finemo/{sample}/report/motif_occurrences.tsv"
    params:
        hits_dir = rules.finemo_call_hits.params.out_dir,
        out_dir  = "finemo/{sample}/report",
        modisco_width = config["modisco_tf"]["window"]
    resources:
        mem_mb  = RESOURCES["finemo_report"]["mem_mb"],
        runtime = RESOURCES["finemo_report"]["runtime"]
    retries:
        RESOURCES["finemo_report"]["retries"]
    conda:
        "../envs/finemo.yaml"
    log:
        "logs/finemo/{sample}_report.log"
    shell:
        "finemo report -r {input.regions} -H {params.hits_dir} -m {input.modisco} -o {params.out_dir} "
        "-W {params.modisco_width} > {log} 2>&1"
