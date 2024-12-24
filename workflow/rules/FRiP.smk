
rule FRiP_stat:
    input:
        peak = "results/callpeak/macs/{peak_type}/{peak}_peaks.{peak_type}Peak",
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

