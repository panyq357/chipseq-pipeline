
ruleorder: macs_with_input > macs_no_input


rule macs_with_input:
    input:
        treatment = lambda w: config["peaks"][w.peak_type][w.peak]["treatment"],
        control = lambda w: config["peaks"][w.peak_type][w.peak]["control"]
    output:
        peak = "results/macs_callpeak/{peak_type}/{peak}_peaks.{peak_type}Peak",
        xls = "results/macs_callpeak/{peak_type}/{peak}_peaks.xls"
    params:
        bam_format = lambda w: config["peaks"][w.peak_type][w.peak]["bam_format"] if "bam_format" in config["peaks"][w.peak_type][w.peak] else config["macs_default_config"]["bam_format"],
        mappable_genome_size = config["mappable_genome_size"],
        name = "{peak}",
        outdir = "results/macs_callpeak/{peak_type}",
        extra_params = lambda w: config["peak_params"][config["peaks"][w.peak_type][w.peak]["extra_params"]] if "extra_params" in config["peaks"][w.peak_type][w.peak] else ""
    resources:
        io = 100
    log:
        "logs/macs_callpeak/{peak_type}/{peak}.log"
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
        peak = "results/macs_callpeak/{peak_type}/{peak}_peaks.{peak_type}Peak",
        xls = "results/macs_callpeak/{peak_type}/{peak}_peaks.xls"
    params:
        bam_format = lambda w: config["peaks"][w.peak_type][w.peak]["bam_format"] if "bam_format" in config["peaks"][w.peak_type][w.peak] else config["macs_default_config"]["bam_format"],
        mappable_genome_size = config["mappable_genome_size"],
        name = "{peak}",
        outdir = "results/macs_callpeak/{peak_type}",
        extra_params = lambda w: config["peak_params"][config["peaks"][w.peak_type][w.peak]["extra_params"]] if "extra_params" in config["peaks"][w.peak_type][w.peak] else ""
    resources:
        io = 100
    log:
        "logs/macs_callpeak/{peak_type}/{peak}.log"
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


rule FRiP_stat:
    input:
        peak = "results/macs_callpeak/{peak_type}/{peak}_peaks.{peak_type}Peak",
        bam = lambda w: config["peaks"][w.peak_type][w.peak]["treatment"]
    output:
        fc = "results/FRiP/{peak}.{peak_type}.fc.rds"
    threads:
        20
    script:
        "../scripts/FRiP_stat.R"


rule FRiP_plot:
    input:
        fc_list = [f"results/FRiP/{peak}.{peak_type}.fc.rds" for peak_type in config["peaks"] for peak in config["peaks"][peak_type]]
    output:
        plot = "results/FRiP.pdf"
    script:
        "../scripts/FRiP_plot.R"

