rule bamCoverage:
    input:
        bam = "{prefix}.bam",
        csi = "{prefix}.bam.csi",
        mappable_genome_size = "results/mappable_genome_size"
    output:
        bw = "{prefix}.bamCoverage.param-{param}.bw"
    params:
        lambda w: config["deeptools_params"]["bamCoverage"][w.param]
    threads:
        2
    resources:
        io = 100
    log:
        "{prefix}.bamCoverage.param-{param}.log"
    shell:
        '''
        bamCoverage \
            --numberOfProcessors {threads} \
            --bam {input.bam} \
            --outFileName {output.bw} \
            --effectiveGenomeSize $(cat {input.mappable_genome_size}) \
            {params} \
            >> {log} 2>&1
        '''


rule computeMatrix:
    input:
        bw = lambda w: config["heatmap_jobs"][w.heatmap]["bw"].values(),
        bed = lambda w: config["heatmap_jobs"][w.heatmap]["bed"].values()
    output:
        temp("results/deeptools/heatmaps/{heatmap}.mat.gz")
    params:
        lambda w: config["deeptools_params"]["computeMatrix"][config["heatmap_jobs"][w.heatmap]["params"]["computeMatrix"]]
    threads:
        20
    log:
        "results/deeptools/heatmaps/{heatmap}.computeMatrix.log"
    resources:
        io = 100
    shell:
        '''
        computeMatrix {params} \
            --numberOfProcessors {threads} \
            --regionsFileName {input.bed} \
            --scoreFileName {input.bw} \
            --outFileName {output} \
            2> {log}
        '''


rule plotHeatmap:
    input:
        "results/deeptools/heatmaps/{heatmap}.mat.gz"
    output:
        "results/deeptools/heatmaps/{heatmap}.heatmap.pdf"
    params:
        sample_labels = lambda w: list(config["heatmap_jobs"][w.heatmap]["bw"].keys()),
        region_labels = lambda w: list(config["heatmap_jobs"][w.heatmap]["bed"].keys()),
        extra_params = lambda w: config["deeptools_params"]["plotHeatmap"][config["heatmap_jobs"][w.heatmap]["params"]["plotHeatmap"]]
    log:
        "results/deeptools/heatmaps/{heatmap}.plotHeatmap.log"
    shell:
        '''
        plotHeatmap {params.extra_params} \
            --matrixFile {input} \
            --outFileName {output} \
            --samplesLabel {params.sample_labels} \
            --regionsLabel {params.region_labels} \
            2> {log}
        '''



rule bigwigAverage:
    input:
        lambda w: config["bigwigAverage_jobs"][w.name]
    output:
        "results/deeptools/bigwigAverage/{name}.bw"
    threads:
        20
    shell:
        "bigwigAverage -b {input} -o {output} -p {threads}"
