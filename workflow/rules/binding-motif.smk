configfile: "config/binding-motif-config.yaml"

def get_peak(peak):
    return f"results/macs_callpeak/{peak}/{peak}_peaks.xls"

rule get_peak_fasta:
    input:
        ref_genome = config["bwa_index_prefix"],
        peak_file = lambda w: get_peak(config["binding-motif-jobs"][w.job]["peak"])
    output:
        out_fasta = "results/get_peak_fasta/{job}.fasta"
    params:
        top = lambda w: config["binding-motif-jobs"][w.job]["top"],
        radius = lambda w: config["binding-motif-jobs"][w.job]["radius"]
    script:
        "../scripts/get_peak_fasta.R"


rule meme:
    input:
        "results/get_peak_fasta/{job}.fasta"
    output:
        directory("results/meme/{job}")
    threads:
        10
    shell:
        "meme-chip -meme-p {threads} --oc {output} {input}"

