
rule get_mappable_genome_size:
    input:
        genome = config["genome"]
    output:
        "results/mappable_genome_size"
    shell:
        '''
        awk '$1 == "total" {{ print $2 }}' <(faCount {input.genome}) > {output}
        '''


rule macs:
    input:
        treatment = lambda w: config["peak_calling_jobs"][w.peak_type][w.peak_name]["treatment"],
        control = branch(
            condition = lambda w: "control" in config["peak_calling_jobs"][w.peak_type][w.peak_name],
            then = lambda w: config["peak_calling_jobs"][w.peak_type][w.peak_name]["control"],
            otherwise = None
        ),
        mappable_genome_size = "results/mappable_genome_size"
    output:
        peak = "results/peak_calling/{peak_type}/{peak_name}_peaks.{peak_type}Peak",
        xls = "results/peak_calling/{peak_type}/{peak_name}_peaks.xls"
    params:
        name = "{peak_name}",
        outdir = "results/peak_calling/{peak_type}",
        extra_params = branch(
            condition = lambda w: "extra_params" in config["peak_calling_jobs"][w.peak_type][w.peak_name],
            then = lambda w: config["macs_params"][config["peak_calling_jobs"][w.peak_type][w.peak_name]["params"]],
            otherwise = config["macs_params"]["default"]
        )
    resources:
        io = 100
    log:
        "results/peak_calling/{peak_type}/{peak_name}_peaks.log"
    run:
        macs_command_lines = [
            "macs2 callpeak",
            "--treatment {input.treatment}",
            "--gsize $(cat {input.mappable_genome_size})",
            "--name {params.name}",
            "--outdir {params.outdir}",
            ">> {log} 2>&1"
        ]
        if input.control is not None:
            macs_command_lines.insert(2, "--control {input.control}")
        if wildcards.peak_type == "broad":
            macs_command_lines.insert(1, "--broad")
        macs_command_lines.insert(2, params.extra_params)
        shell(" ".join(macs_command_lines))


rule FRiP_stat:
    input:
        peak = "results/peak_calling/{peak_type}/{peak_name}_peaks.{peak_type}Peak",
        treatment = lambda w: config["peak_calling_jobs"][w.peak_type][w.peak_name]["treatment"],
    output:
        fc = "results/peak_calling/{peak_type}/{peak_name}_peaks.FRiP.rds"
    wildcard_constraints:
        treatment = "[^.]+",
        format = "[^.]+"
    threads:
        20
    script:
        "../scripts/FRiP_stat.R"


rule FRiP_plot:
    input:
        fc_list = [f"results/peak_calling/{peak_type}/{peak_name}_peaks.FRiP.rds" for peak_type in config["peak_calling_jobs"] for peak_name in config["peak_calling_jobs"][peak_type]]
    output:
        plot = "results/FRiP.pdf"
    script:
        "../scripts/FRiP_plot.R"
