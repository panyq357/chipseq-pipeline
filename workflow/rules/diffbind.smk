import pandas as pd

configfile: "config/diffbind-config.yaml"


def get_input_files_from_sample_sheet(sample_sheet_path):
    sample_sheet = pd.read_csv(sample_sheet_path)
    input_files = sample_sheet["bamReads"].to_list() + sample_sheet["Peaks"].to_list()
    if "bamControl" in sample_sheet.columns:
        input_files += sample_sheet["bamControl"].to_list()
    return input_files


rule diffbind_counts:
    input:
        lambda w: get_input_files_from_sample_sheet(config["diffbind-sample-sheets"][w.sheet]),
        sample_sheet = lambda w: config["diffbind-sample-sheets"][w.sheet]
    output:
        counts_rds = "results/diffbind_counts/{sheet}_counts.rds",
        peakset_csv = "results/diffbind_counts/{sheet}_peakset.csv"
    threads:
        6  # Less cores, less memory consumption. For 6 cores, it needs about 30G memory.
    script:
        "../scripts/diffbind_counts.R"


rule diffbind_analyze:
    input:
        counts_rds = "results/diffbind_counts/{sheet}_counts.rds"
    output:
        ma_plot = "results/diffbind_analyze/{sheet}/{sheet}.ma_plot.svg",
        volcano = "results/diffbind_analyze/{sheet}/{sheet}.volcano.svg",
        report_rds = "results/diffbind_analyze/{sheet}/{sheet}.report.rds",
        report_csv = "results/diffbind_analyze/{sheet}/{sheet}.report.csv",
        fdr_peak = "results/diffbind_analyze/{sheet}/{sheet}.fdr_peak.bed",
        fdr_peak_up = "results/diffbind_analyze/{sheet}/{sheet}.fdr_peak_up.bed",
        fdr_peak_down = "results/diffbind_analyze/{sheet}/{sheet}.fdr_peak_down.bed"
    params:
        fdr_th = 1,
        p_th = 0.05
    script:
        "../scripts/diffbind_analyze.R"


rule diff_peak_matrix:
    input:
        fdr_peak_up = "results/diffbind_analyze/{sheet}/{sheet}.fdr_peak_up.bed",
        fdr_peak_down = "results/diffbind_analyze/{sheet}/{sheet}.fdr_peak_down.bed",
        bw = lambda w: config["diff-peak-heatmap-bw"][w.sheet]
    output:
        temp("results/diff_peak_heatmap/{sheet}.diff_peak_matrix.mat.gz")
    threads:
        20
    shell:
        '''
        computeMatrix reference-point \
            -p {threads} \
            --referencePoint center \
            --upstream 1000 \
            --downstream 1000 \
            --regionsFileName {input.fdr_peak_up} {input.fdr_peak_down} \
            --scoreFileName {input.bw} \
            --outFileName {output}
        '''

rule diff_peak_heatmap:
    input:
        "results/diff_peak_heatmap/{sheet}.diff_peak_matrix.mat.gz"
    output:
        "results/diff_peak_heatmap/{sheet}.diff_peak_heatmap.png"
    shell:
        '''
        plotHeatmap \
            --legendLocation "upper-right" \
            --zMax 3 \
            --xAxisLabel "" \
            --regionsLabel "Up" "Down" \
            -m {input} \
            -out {output}
        '''


