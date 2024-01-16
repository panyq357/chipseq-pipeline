configfile: "config/reads2bwa-config.yaml"


# paired end is prefered over single end.
ruleorder: bwa_mapping_pe > bwa_mapping_se > fastp_pe > fastp_se > cat_pe > cat_se


# Rule cat_pe and cat_se is for unifying format of raw FASTQ files,
# All FASTQ will be gzipped and concated before sending to fastp.
rule cat_pe:
    input:
        r1 = lambda w: config["samples"][w.sample_id]["R1"],
        r2 = lambda w: config["samples"][w.sample_id]["R2"]
    output:
        r1_cat = temp("resources/cat/{sample_id}.R1.fastq.gz"),
        r2_cat = temp("resources/cat/{sample_id}.R2.fastq.gz")
    log:
        "logs/cat/{sample_id}.log"
    run:
        if type(input.r1) is str:  # Only one file.
            if input.r1[-3:] == ".gz":
                shell("ln -s $(realpath {input.r1}) {output.r1_cat} 2> {log}")
            else:
                shell("gzip -c {input.r1} > {output.r1_cat} 2> {log}")
        else:  # Multiple files.
            for f in input.r1:
                if f[-3:] == ".gz":
                    shell("cat {f} >> {output.r1_cat} 2> {log}")
                else:
                    shell("gzip -c {f} >> {output.r1_cat} 2> {log}")
        if type(input.r2) is str:
            if input.r2[-3:] == ".gz":
                shell("ln -s $(realpath {input.r2}) {output.r2_cat} 2> {log}")
            else:
                shell("gzip -c {input.r2} > {output.r2_cat} 2> {log}")
        else:
            for f in input.r2:
                if f[-3:] == ".gz":
                    shell("cat {f} >> {output.r2_cat} 2> {log}")
                else:
                    shell("gzip -c {f} >> {output.r2_cat} 2> {log}")


rule cat_se:
    input:
        r1 = lambda w: config["samples"][w.sample_id]["R1"],
    output:
        r1_cat = temp("resources/cat/{sample_id}.R1.fastq.gz")
    log:
        "logs/cat/{sample_id}.log"
    run:
        if type(input.r1) is str:  # Only one file.
            if input.r1[-3:] == ".gz":
                shell("ln -s $(realpath {input.r1}) {output.r1_cat} 2> {log}")
            else:
                shell("gzip -c {input.r1} > {output.r1_cat} 2> {log}")
        else:  # Multiple files.
            for f in input.r1:
                if f[-3:] == ".gz":
                    shell("cat {f} >> {output.r1_cat} 2> {log}")
                else:
                    shell("gzip -c {f} >> {output.r1_cat} 2> {log}")


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
        config["threads"]["fastp_pe"]
    priority:
        60
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
        config["threads"]["fastp_se"]
    priority:
        60
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
        config["threads"]["bwa_mapping_pe"]
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
        config["threads"]["bwa_mapping_se"]
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

