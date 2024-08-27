
# paired end is prefered over single end.
ruleorder: bwa_mapping_pe > bwa_mapping_se


rule bwa_mapping_pe:
    input:
        r1_fastp = rules.fastp_pe.output.r1_fastp,
        r2_fastp = rules.fastp_pe.output.r2_fastp,
        bwa_index = [
            f"{config['bwa_index_prefix']}.{suffix}" for suffix in ["amb", "ann", "bwt", "pac", "sa"]
        ]
    output:
        bam = "results/bwa_mapping/{sample_id}.bam",
        bai = "results/bwa_mapping/{sample_id}.bam.bai"
    params:
        bwa_index_prefix = config['bwa_index_prefix'],
        samtools_sort_threads = 4  # Less threads, less memory.
    log:
        "logs/bwa_mapping/{sample_id}.log"
    threads:
        20
    priority:
        70
    shell:
        '''
        bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.sample_id}\\tSM:{wildcards.sample_id}' \
            {params.bwa_index_prefix} {input.r1_fastp} {input.r2_fastp} 2>> {log} \
        | samtools fixmate -u -m - - 2>> {log} \
        | samtools sort -u -@{params.samtools_sort_threads} 2>> {log} \
        | samtools markdup -u -@{threads} - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        '''


rule bwa_mapping_se:
    input:
        r1_fastp = rules.fastp_se.output.r1_fastp,
        bwa_index = [
            f"{config['bwa_index_prefix']}.{suffix}" for suffix in ["amb", "ann", "bwt", "pac", "sa"]
        ]
    output:
        bam = "results/bwa_mapping/{sample_id}.bam",
        bai = "results/bwa_mapping/{sample_id}.bam.bai"
    params:
        bwa_index_prefix = config['bwa_index_prefix']
    log:
        "logs/bwa_mapping/{sample_id}.log"
    threads:
        20
    priority:
        70
    shell:
        '''
        bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.sample_id}\\tSM:{wildcards.sample_id}' \
            {params.bwa_index_prefix} {input.r1_fastp} 2>> {log} \
        | samtools fixmate -u -m - - 2>> {log} \
        | samtools sort -u -@{threads} 2>> {log} \
        | samtools markdup -u -@{threads} - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        '''


rule flagstats:
    input:
        bam = "results/bwa_mapping/{sample_id}.bam",
        bai = "results/bwa_mapping/{sample_id}.bam.bai"
    output:
        txt = "results/bwa_mapping/{sample_id}.flagstats.txt"
    log:
        "logs/flagstats/{sample_id}.log"
    priority:
        80
    shell:
        '''
        samtools flagstats {input.bam} > {output.txt} 2> {log}
        '''

