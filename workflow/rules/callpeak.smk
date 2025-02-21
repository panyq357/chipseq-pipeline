
ruleorder: macs_with_input > macs_no_input


rule macs_with_input:
    input:
        treatment = lambda w: config["peaks"][w.peak_type][w.peak]["treatment"],
        control = lambda w: config["peaks"][w.peak_type][w.peak]["control"]
    output:
        peak = "results/callpeak/macs/{peak_type}/{peak}_peaks.{peak_type}Peak",
        xls = "results/callpeak/macs/{peak_type}/{peak}_peaks.xls"
    params:
        bam_format = lambda w: config["peaks"][w.peak_type][w.peak]["bam_format"] if "bam_format" in config["peaks"][w.peak_type][w.peak] else config["macs_default_config"]["bam_format"],
        mappable_genome_size = config["mappable_genome_size"],
        name = "{peak}",
        outdir = "results/callpeak/macs/{peak_type}",
        extra_params = lambda w: config["peak_params"][config["peaks"][w.peak_type][w.peak]["extra_params"]] if "extra_params" in config["peaks"][w.peak_type][w.peak] else ""
    resources:
        io = 100
    log:
        "logs/callpeak/macs/{peak_type}/{peak}.log"
    run:
        macs_command_lines = [
            "macs2 callpeak",
            "--treatment {input.treatment} ",
            "--control {input.control} ",
            "--format {params.bam_format}",
            "--gsize {params.mappable_genome_size}",
            "--name {params.name}",
            "--outdir {params.outdir}",
            ">> {log} 2>&1"
        ]
        if wildcards.peak_type == "broad":
            macs_command_lines.insert(1, "--broad")
        macs_command_lines.insert(2, params.extra_params)
        shell(" ".join(macs_command_lines))


rule macs_no_input:
    input:
        treatment = lambda w: config["peaks"][w.peak_type][w.peak]["treatment"],
    output:
        peak = "results/callpeak/macs/{peak_type}/{peak}_peaks.{peak_type}Peak",
        xls = "results/callpeak/macs/{peak_type}/{peak}_peaks.xls"
    params:
        bam_format = lambda w: config["peaks"][w.peak_type][w.peak]["bam_format"] if "bam_format" in config["peaks"][w.peak_type][w.peak] else config["macs_default_config"]["bam_format"],
        mappable_genome_size = config["mappable_genome_size"],
        name = "{peak}",
        outdir = "results/callpeak/macs/{peak_type}",
        extra_params = lambda w: config["peak_params"][config["peaks"][w.peak_type][w.peak]["extra_params"]] if "extra_params" in config["peaks"][w.peak_type][w.peak] else ""
    resources:
        io = 100
    log:
        "logs/callpeak/macs/{peak_type}/{peak}.log"
    run:
        macs_command_lines = [
            "macs2 callpeak",
            "--treatment {input.treatment} ",
            "--format {params.bam_format}",
            "--gsize {params.mappable_genome_size}",
            "--name {params.name}",
            "--outdir {params.outdir}",
            "--nolambda",  # better background control for peak calling without control sample.
            ">> {log} 2>&1"
        ]
        if wildcards.peak_type == "broad":
            macs_command_lines.insert(1, "--broad")
        macs_command_lines.insert(2, params.extra_params)
        shell(" ".join(macs_command_lines))
