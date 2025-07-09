'''
Unify format of raw FASTQ files,
All FASTQ will be gzipped and concated before sending to fastp.
'''


rule concat_fastq_paired:
    input:
        lambda w: config["read_mapping_jobs"]["paired"][w.sample_id][w.end],
    output:
        temp("resources/concat_fastq/paired/{sample_id}.{end}.fastq.gz"),
    log:
        "resources/concat_fastq/paired/{sample_id}.{end}.log"
    resources:
        io = 50
    run:
        # If only one gzipped file, soft link it to save IO.
        if len(input) == 1 and input[0].lower().endswith(".gz"):
            shell("ln -s $(realpath {input[0]}) {output} 2> {log}")
            return()

        for f in input:
            if f.endswith(".gz"):
                shell("cat {f} >> {output} 2>> {log}")
            elif f.endswith(".bz2"):
                shell("bzip2 -dc {f} 2>> {log} | gzip -c >> {output} 2>> {log}")
            elif f.endswith(".fasta") or f.endswith(".fa"):
                shell("gzip -c {f} >> {output} 2>> {log}")


rule concat_fastq_single:
    input:
        lambda w: config["read_mapping_jobs"]["single"][w.sample_id],
    output:
        temp("resources/concat_fastq/single/{sample_id}.single.fastq.gz"),
    log:
        "resources/concat_fastq/single/{sample_id}.single.log"
    resources:
        io = 50
    run:
        # If only one gzipped file, soft link it to save IO.
        if len(input) == 1 and input[0].lower().endswith(".gz"):
            shell("ln -s $(realpath {input[0]}) {output} 2> {log}")
            return()

        for f in input:
            if f.endswith(".gz"):
                shell("cat {f} >> {output} 2>> {log}")
            elif f.endswith(".bz2"):
                shell("bzip2 -dc {f} 2>> {log} | gzip -c >> {output} 2>> {log}")
            elif f.endswith(".fasta") or f.endswith(".fa"):
                shell("gzip -c {f} >> {output} 2>> {log}")


'''
Use fastp to filter FASTQ, and generate report html.
'''


rule fastp_paired:
    input:
        r1 = "resources/concat_fastq/paired/{sample_id}.r1.fastq.gz",
        r2 = "resources/concat_fastq/paired/{sample_id}.r2.fastq.gz",
    output:
        r1 = temp("resources/fastp/paired/{sample_id}.fastp.r1.fastq.gz"),
        r2 = temp("resources/fastp/paired/{sample_id}.fastp.r2.fastq.gz"),
        json = "results/fastp/paired/{sample_id}.fastp.json",
        html = "results/fastp/paired/{sample_id}.fastp.html"
    log:
        "results/fastp/paired/{sample_id}.fastp.log"
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


rule fastp_single:
    input:
        se = "resources/concat_fastq/single/{sample_id}.single.fastq.gz",
    output:
        se = temp("resources/fastp/single/{sample_id}.fastp.single.fastq.gz"),
        json = "results/fastp/single/{sample_id}.fastp.single.json",
        html = "results/fastp/single/{sample_id}.fastp.single.html"
    log:
        "results/fastp/single/{sample_id}.fastp.log"
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
            --in1 {input.se} \
            --out1 {output.se} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample_id} \
            --verbose 2>> {log}
        '''


'''
Use BWA MEM to map reads to reference.
'''


rule bwa_mem_paired:
    input:
        r1 = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_paired.input.r1,
            otherwise = rules.fastp_paired.output.r1
        ),
        r2 = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_paired.input.r2,
            otherwise = rules.fastp_paired.output.r2
        ),
        bwa_index = multiext(config["genome"],".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = "results/read_mapping/paired/{sample_id}.bam"
    params:
        index_prefix = config["genome"],
        read_group = "@RG\\tID:{sample_id}\\tSM:{sample_id}",
        sort_threads = config["samtools"]["sort"]["threads"],
        sort_mem_per_thread = config["samtools"]["sort"]["mem_per_thread"]
    log:
        "results/read_mapping/paired/{sample_id}.log"
    wildcard_constraints:
        sample_id ="[^.]+"
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
        se = branch(
            condition = config["fastp"]["skip"],
            then = rules.fastp_single.input.se,
            otherwise = rules.fastp_single.output.se
        ),
        bwa_index = multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = "results/read_mapping/single/{sample_id}.bam"
    params:
        index_prefix = config["genome"],
        read_group = "@RG\\tID:{sample_id}\\tSM:{sample_id}",
        sort_threads = config["samtools"]["sort"]["threads"],
        sort_mem_per_thread = config["samtools"]["sort"]["mem_per_thread"]
    log:
        "results/read_mapping/single/{sample_id}.log"
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
        "{prefix}.flagstats.txt"
    threads:
        4
    shell:
        "samtools flagstats -@ {threads} {input.bam} > {output}"


rule collect_all_flagstats:
    input:
        [f"results/read_mapping/{end}/{name}.flagstats.txt" for end, names in config["read_mapping_jobs"].items() for name in names]
    output:
        "results/all_flagstats.tsv"
    script:
        "../scripts/collect_all_flagstats.R"


rule bam_filter:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.filter-{filter}.bam",
    params:
        lambda w: config["bam_filters"][w.filter]
    resources:
        io = 50
    threads:
        4
    shell:
        "samtools view -@{threads} -b -h {params} {input} > {output}"

