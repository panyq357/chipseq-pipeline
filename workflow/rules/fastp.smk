
'''
Use fastp to filter FASTQ, and generate report html.
'''

rule fastp_pe:
    input:
        r1 = "resources/cat/{sample_id}.R1.fastq.gz",
        r2 = "resources/cat/{sample_id}.R2.fastq.gz",
    output:
        r1 = temp("resources/fastp/{sample_id}.fastp.R1.fastq.gz"),
        r2 = temp("resources/fastp/{sample_id}.fastp.R2.fastq.gz"),
        json = "results/fastp/{sample_id}.fastp.json",
        html = "results/fastp/{sample_id}.fastp.html"
    log:
        "results/fastp/{sample_id}.fastp.log"
    threads:
        4
    priority:
        60
    resources:
        io = 100
    shell:
        '''
        fastp \
            --thread {threads} \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''


rule fastp_se:
    input:
        r1 = "resources/cat/{sample_id}.R1.fastq.gz",
    output:
        se = temp("resources/fastp/{sample_id}.fastp.fastq.gz"),
        json = "results/fastp/{sample_id}.fastp.json",
        html = "results/fastp/{sample_id}.fastp.html"
    log:
        "results/fastp/{sample_id}.fastp.log"
    threads:
        4
    priority:
        60
    resources:
        io = 100
    shell:
        '''
        fastp \
            --thread {threads} \
            --in1 {input.r1} \
            --out1 {output.se} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''

