
'''
Use fastp to filter FASTQ, and generate report html.
'''

ruleorder: fastp_pe > fastp_se

rule fastp_pe:
    input:
        r1 = rules.cat_pe.output.r1_cat,
        r2 = rules.cat_pe.output.r2_cat
    output:
        r1_fastp = "resources/fastp/{sample_id}.fastp.R1.fastq.gz",
        r2_fastp = "resources/fastp/{sample_id}.fastp.R2.fastq.gz",
        json = "results/fastp/{sample_id}.json",
        html = "results/fastp/{sample_id}.html"
    log:
        "logs/fastp/{sample_id}.log"
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
            --out1 {output.r1_fastp} \
            --out2 {output.r2_fastp} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''


rule fastp_se:
    input:
        r1 = rules.cat_se.output.r1_cat
    output:
        r1_fastp = "resources/fastp/{sample_id}.fastp.R1.fastq.gz",
        json = "results/fastp/{sample_id}.json",
        html = "results/fastp/{sample_id}.html"
    log:
        "logs/fastp/{sample_id}.log"
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
            --out1 {output.r1_fastp} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''

