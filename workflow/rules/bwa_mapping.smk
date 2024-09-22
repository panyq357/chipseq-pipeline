
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
        sort_threads = config["samtools"]["sort_threads"],
        sort_mem_per_thread = config["samtools"]["sort_mem_per_thread"]
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
        | samtools sort -u -@{params.sort_threads} -m{params.sort_mem_per_thread} 2>> {log} \
        | samtools markdup -u -@{threads} - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        samtools index -@{threads} {output.bam} 2>> {log}
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
        bwa_index_prefix = config['bwa_index_prefix'],
        sort_threads = config["samtools"]["sort_threads"],
        sort_mem_per_thread = config["samtools"]["sort_mem_per_thread"]
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
        | samtools sort -u -@{params.sort_threads} -m{params.sort_mem_per_thread} 2>> {log} \
        | samtools markdup -u -@{threads} - - 2>> {log} \
        | samtools view -b -h -@{threads} -o {output.bam} 2>> {log}
        samtools index -@{threads} {output.bam} 2>> {log}
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
    threads:
        4
    shell:
        '''
        samtools flagstats -@{input.bam} > {output.txt} 2> {log}
        '''


rule bam_filter:
    input:
        bam = "results/bwa_mapping/{sample_id}.bam",
        bai = "results/bwa_mapping/{sample_id}.bam.bai"
    output:
        bam = "results/bam_filter/{sample_id}.filter.bam",
        bai = "results/bam_filter/{sample_id}.filter.bam.bai"
    params:
        filter_args = config["samtools"]["filter_args"]
    resources:
        io = 50
    threads:
        4
    shell:
        '''
        samtools view -@{threads} -b -h {params.filter_args} {input.bam} > {output.bam}
        samtools index -@{threads} {output.bam}
        '''

rule bam_merge:
    input:
        bams = lambda w: config["bam_merge"][w.name]["bams"]
    output:
        bam = "results/bam_merge/{name}.bam",
        bai = "results/bam_merge/{name}.bam.bai"
    resources:
        io = 50
    threads:
        4
    shell:
        '''
        samtools merge -@{threads} {output.bam} {input.bams}
        samtools index -@{threads} {output.bam}
        '''
