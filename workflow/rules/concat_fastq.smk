
'''
Unify format of raw FASTQ files,
All FASTQ will be gzipped and concated before sending to fastp.
'''

rule concat_fastq:
    input:
        lambda w: config["samples"][w.sample_id][w.end],
    output:
        temp("resources/cat/{sample_id}.{end}.fastq.gz"),
    log:
        "resources/cat/{sample_id}.{end}.log"
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

