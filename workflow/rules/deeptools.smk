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


rule bigwigAverage:
    input:
        lambda w: config["bigwigAverage_jobs"][w.name]
    output:
        "results/deeptools/bigwigAverage/{name}.bw"
    threads:
        20
    shell:
        "bigwigAverage -b {input} -o {output} -p {threads}"


rule computeMatrix:
    input:
        bw = lambda w: config["computeMatrix_jobs"][w.name]["bw"].values(),
        bed = lambda w: config["computeMatrix_jobs"][w.name]["bed"].values()
    output:
        "results/deeptools/plots/{name}.mat.gz"
    params:
        lambda w: config["deeptools_params"]["computeMatrix"][config["computeMatrix_jobs"][w.name]["param"]]
    threads:
        20
    log:
        "results/deeptools/plots/{name}.log"
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
        "results/deeptools/plots/{name}.mat.gz"
    output:
        "results/deeptools/plots/{name}.heatmap-{param}.pdf"
    params:
        sample_labels = lambda w: list(config["computeMatrix_jobs"][w.name]["bw"].keys()),
        region_labels = lambda w: list(config["computeMatrix_jobs"][w.name]["bed"].keys()),
        extra_params = lambda w: config["deeptools_params"]["plotHeatmap"][w.param]
    log:
        "results/deeptools/plots/{name}.heatmap-{param}.log"
    wildcard_constraints:
        name = "[^.]+"
    shell:
        '''
        plotHeatmap {params.extra_params} \
            --matrixFile {input} \
            --outFileName {output} \
            --samplesLabel {params.sample_labels} \
            --regionsLabel {params.region_labels} \
            {params.extra_params}
            2> {log}
        '''




rule plotProfile:
    input:
        "results/deeptools/plots/{name}.mat.gz"
    output:
        "results/deeptools/plots/{name}.profile-{param}.pdf"
    params:
        sample_labels = lambda w: list(config["computeMatrix_jobs"][w.name]["bw"].keys()),
        region_labels = lambda w: list(config["computeMatrix_jobs"][w.name]["bed"].keys()),
        extra_params = lambda w: config["deeptools_params"]["plotProfile"][w.param]
    log:
        "results/deeptools/plots/{name}.profile-{param}.log"
    shell:
        '''
        plotProfile \
            --matrixFile {input} \
            --outFileName {output} \
            --samplesLabel {params.sample_labels} \
            --regionsLabel {params.region_labels} \
            {params.extra_params}
            2> {log}
        '''
