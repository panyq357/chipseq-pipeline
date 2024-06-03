configfile: "config/deeptools-config.yaml"


rule bamCoverage:
    input:
        bam = "results/bwa_mapping/{sample}.bam"
    output:
        bw = "results/bamCoverage/{sample}.bw"
    params:
        mappable_genome_size = config["mappable_genome_size"]
    threads:
        2
    resources:
        disk_io = 10
    log:
        "logs/bamCoverage/{sample}.log"
    shell:
        '''
        bamCoverage \
            -p {threads} \
            --bam {input.bam} \
            --effectiveGenomeSize {params.mappable_genome_size} \
            --normalizeUsing CPM \
            --outFileName {output.bw} >> {log} 2>&1
        '''


def get_bam(sample):
    return f"results/bwa_mapping/{sample}.bam"


rule bamCompare:
    input:
        b1 = lambda w: get_bam(config["peaks"][w.peak_with_input]["ip"]),
        b2 = lambda w: get_bam(config["peaks"][w.peak_with_input]["input"])
    output:
        "results/bamCompare/{peak_with_input}-compare.bw"
    params:
        mappable_genome_size = config["mappable_genome_size"]
    threads:
        2
    resources:
        disk_io = 10
    log:
        "logs/bamCompare/{peak_with_input}.log"
    shell:
        '''
        bamCompare \
            -p {threads} \
            -b1 {input.b1} \
            -b2 {input.b2} \
            --effectiveGenomeSize {params.mappable_genome_size} \
            --scaleFactorsMethod None \
            --normalizeUsing CPM \
            --outFileName {output} >> {log} 2>&1
        '''


rule tss_to_tes_matrix:
    input:
        bw = lambda w: config["tss-to-tes-jobs"][w.job],
        bed = config["transcript_bed"]
    output:
        temp("results/tss_to_tes_heatmap/{job}.tss_to_tes.matrix.mat.gz")
    threads:
        20
    shell:
        '''
        computeMatrix scale-regions \
            -p {threads} \
            --regionBodyLength 5000 \
            --upstream 2000 \
            --downstream 2000 \
            --regionsFileName {input.bed} \
            --scoreFileName {input.bw} \
            --outFileName {output}
        '''


rule tss_to_tes_heatmap:
    input:
        "results/tss_to_tes_heatmap/{job}.tss_to_tes.matrix.mat.gz"
    output:
        "results/tss_to_tes_heatmap/{job}.tss_to_tes.heatmap.png"
    shell:
        '''
        plotHeatmap \
            --legendLocation none \
            --zMax 7 \
            --xAxisLabel "" \
            --regionsLabel "" \
            -m {input} \
            -out {output}
        '''


rule tss_to_tes_profile:
    input:
        "results/tss_to_tes_heatmap/{job}.tss_to_tes.matrix.mat.gz"
    output:
        "results/tss_to_tes_heatmap/{job}.tss_to_tes.profile.png"
    shell:
        '''
        plotProfile -m {input} \
            -out {output} \
            --perGroup \
            --plotTitle ""
        '''


def get_peak_bw(peak):

    def get_sample_bw(sample):
        return f"results/bamCoverage/{sample}.bw"

    ip_sample = config["peaks"][peak]["ip"]
    if "input" in config["peaks"][peak]:
        input_sample = config["peaks"][peak]["input"]
        return [get_sample_bw(ip_sample), get_sample_bw(input_sample)]
    else:
        return get_sample_bw(ip_sample)


rule peak_heatmap:
    input:
        bw = lambda w: get_peak_bw(w.peak)
    output:
        mat = temp("results/peak_heatmap/{peak}.matrix.mat.gz"),
        png = "results/peak_heatmap/{peak}.heatmap.png"
    params:
        peak_bed = lambda w: get_peak_bed(w.peak)
    threads:
        20
    log:
        "logs/peak_heatmap/{peak}.log"
    shell:
        '''
        computeMatrix reference-point \
            -p {threads} \
            --referencePoint center \
            --upstream 1000 \
            --downstream 1000 \
            --regionsFileName {params.peak_bed} \
            --scoreFileName \
                {input.bw} \
            --outFileName {output.mat} >> {log} 2>&1
        plotHeatmap \
            -m {output.mat} \
            --legendLocation none \
            -out {output.png} >> {log} 2>&1
        '''


def get_peak_bed(peak):
    return f"results/macs_callpeak/{peak}/{peak}_peaks." + config["peaks"][peak]["type"] + "Peak"


rule FRiP_stat:
    input:
        peak = "results/macs_callpeak/{peak}/{peak}_peaks.xls",
        bam = lambda w: get_bam(config["peaks"][w.peak]["ip"])
    params:
        peak_bed = lambda w: get_peak_bed(w.peak)
    output:
        fc = "results/FRiP/{peak}.fc.rds"
    threads:
        20
    script:
        "../scripts/FRiP_stat.R"


rule FRiP_plot:
    input:
        fc_list = expand("results/FRiP/{peak}.fc.rds", peak=config["peaks"])
    output:
        plot = "results/FRiP.pdf"
    script:
        "../scripts/FRiP_plot.R"

