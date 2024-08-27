
'''
Unify format of raw FASTQ files,
All FASTQ will be gzipped and concated before sending to fastp.
'''

ruleorder: cat_pe > cat_se

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


