
rule slice_fa:
    input:
        peak = "{prefix}.xls",
        genome = config["genome"]
    output:
        fa = "{prefix}.sort-by-{by}-{from}-{to}-{mode}-radius-{radius}.fa"
    script:
        "../scripts/slice_fa.R"


rule meme:
    input:
        fa = "{prefix}.fa",
        db = config["meme"]["database"]
    output:
        directory("{prefix}.meme-chip")
    threads:
        10
    run:
        cmd = "meme-chip -meme-p {threads} --oc {output} {input}"
        for x in input.db:
            cmd = cmd + " -db " + x
        print(cmd)
        shell(cmd)


