

rule get_peak_fasta:
    input:
        ref_genome = config["bwa_index_prefix"],
        peak_xls = lambda w: config["binding-motif-jobs"][w.job]["peak_xls"]
    output:
        out_fasta = "results/get_peak_fasta/{job}.fasta"
    params:
        lambda w: config["binding-motif-jobs"][w.job]["params"],
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

