rule slice_fa:
    input:
        peak = "results/peak_calling/{peak_type}/{peak_name}_peaks.{peak_type}Peak",
        genome = config["genome"]
    output:
        fa = "results/peak_calling/{peak_type}/{peak_name}_peaks.sort-by-{by}-{from}-{to}-{mode}-radius-{radius}.fa"
    script:
        "../scripts/slice_fa.R"


rule meme:
    input:
        fa = "{prefix}.fa",
        db = config["meme_db"]
    output:
        directory("{prefix}.meme-chip")
    threads:
        10
    run:
        cmd = "meme-chip -meme-p {threads} --oc {output} {input}"
        for x in input.db:
            cmd = cmd + " -db " + x
        shell(cmd)
