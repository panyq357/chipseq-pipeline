
rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        "samtools index {input}"


rule flagstats:
    input:
        bam = "{prefix}.bam",
        bai = "{prefix}.bam.bai"
    output:
        txt = "{prefix}.flagstats.txt"
    threads:
        4
    shell:
        "samtools flagstats -@ {threads} {input.bam} > {output.txt}"


ruleorder: bam_filter > bwa_mapping_pe > bwa_mapping_se


rule bam_filter:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.filter.bam",
    params:
        config["samtools"]["filter_params"]
    resources:
        io = 50
    threads:
        4
    shell:
        "samtools view -@{threads} -b -h {params} {input} > {output}"
