
'''
Use BWA MEM to map reads to reference.
'''

rule bwa_mem_paired:
    input:
        r1 = "resources/cat/{sample_id}.R1.fastq.gz" if config["skip_fastp"] else "resources/fastp/{sample_id}.fastp.R1.fastq.gz",
        r2 = "resources/cat/{sample_id}.R2.fastq.gz" if config["skip_fastp"] else "resources/fastp/{sample_id}.fastp.R2.fastq.gz",
        bwa_index = multiext(config["genome"],".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = "results/read_mapping/bwa_mem/{sample_id}.paired.bam"
    params:
        index_prefix = config["genome"],
        read_group = "@RG\\tID:{sample_id}\\tSM:{sample_id}",
        sort_threads = config["samtools"]["sort"]["threads"],
        sort_mem_per_thread = config["samtools"]["sort"]["mem_per_thread"]
    log:
        "results/read_mapping/bwa_mem/{sample_id}.paired.log"
    threads:
        20
    priority:
        70
    shell:
        '''
        bwa mem -M -t {threads} -R '{params.read_group}' {params.index_prefix} {input.r1} {input.r2} 2>> {log} \
        | samtools fixmate -u -m - - 2>> {log} \
        | samtools sort -u -@ {params.sort_threads} -m {params.sort_mem_per_thread} 2>> {log} \
        | samtools markdup -u - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        '''


rule bwa_mem_single:
    input:
        se = "resources/cat/{sample_id}.SE.fastq.gz" if config["skip_fastp"] else "resources/fastp/{sample_id}.fastp.SE.fastq.gz",
        bwa_index = multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = "results/read_mapping/bwa_mem/{sample_id}.single.bam"
    params:
        index_prefix = config["genome"],
        read_group = "@RG\\tID:{sample_id}\\tSM:{sample_id}",
        sort_threads = config["samtools"]["sort"]["threads"],
        sort_mem_per_thread = config["samtools"]["sort"]["mem_per_thread"]
    log:
        "results/read_mapping/bwa_mem/{sample_id}.single.log"
    threads:
        20
    priority:
        70
    shell:
        '''
        bwa mem -M -t {threads} -R {params.read_group} {params.index_prefix} {input.se} 2>> {log} \
        | samtools fixmate -u -m - - 2>> {log} \
        | samtools sort -u -@ {params.sort_threads} -m {params.sort_mem_per_thread} 2>> {log} \
        | samtools markdup -u - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        '''

