

ruleorder: bamCompare > bamCoverage > bigwigAverage


rule bamCoverage:
    input:
        bam = lambda w: config["bws"][w.sample]["b1"],
        bai = lambda w: config["bws"][w.sample]["b1"] + ".bai"
    output:
        bw = "results/deeptools/bw/{sample}.bw"
    params:
        mappable_genome_size = config["mappable_genome_size"],
        bamCoverage_args = config["deeptools"]["bamCoverage"]
    threads:
        2
    resources:
        io = 100
    log:
        "logs/deeptools/bw/{sample}.log"
    shell:
        '''
        bamCoverage \
            --numberOfProcessors {threads} \
            --bam {input.bam} \
            --outFileName {output.bw} \
            --effectiveGenomeSize {params.mappable_genome_size} \
            {params.bamCoverage_args} \
            >> {log} 2>&1
        '''


rule bamCompare:
    input:
        b1 = lambda w: config["bws"][w.sample]["b1"],
        b2 = lambda w: config["bws"][w.sample]["b2"],
        b1_bai = lambda w: config["bws"][w.sample]["b1"] + ".bai",
        b2_bai = lambda w: config["bws"][w.sample]["b2"] + ".bai"
    output:
        bw = "results/deeptools/bw/{sample}.bw"
    params:
        mappable_genome_size = config["mappable_genome_size"],
        bamCompare_args = config["deeptools"]["bamCompare"]
    threads:
        2
    resources:
        io = 100
    log:
        "logs/deeptools/bw/{sample}.log"
    shell:
        '''
        bamCompare \
            --numberOfProcessors {threads} \
            --bamfile1 {input.b1} \
            --bamfile2 {input.b2} \
            --outFileName {output.bw} \
            --effectiveGenomeSize {params.mappable_genome_size} \
            {params.bamCompare_args} 2> {log}
        '''

rule bigwigAverage:
    input:
        bam_list = lambda w: config["bws"][w.sample]["average"],
        bai_list = lambda w: [x + ".bai" for x in config["bws"][w.sample]["average"]]
    output:
        bw = "results/deeptools/bw/{sample}.bw"
    threads:
        2
    resources:
        io = 100
    log:
        "logs/deeptools/bw/{sample}.log"
    shell:
        '''
        bigwigAverage -b {input.bam_list} -o {output.bw} 2> {log}
        '''


rule computeMatrix:
    input:
        bw = lambda w: config["deeptools_heatmaps"][w.heatmap]["bw"],
        bed = lambda w: config["deeptools_heatmaps"][w.heatmap]["bed"]
    output:
        temp("results/deeptools/heatmaps/{heatmap}.mat.gz")
    params:
        lambda w: config["deeptools_params"][config["deeptools_heatmaps"][w.heatmap]["params"]]["matrix"]
    threads:
        20
    log:
        "logs/deeptools/heatmaps/{heatmap}.computeMatrix.log"
    shell:
        '''
        computeMatrix {params} \
            --numberOfProcessors {threads} \
            --regionsFileName {input.bed} \
            --scoreFileName {input.bw} \
            --outFileName {output} \
            > {log}
        '''


rule plotHeatmap:
    input:
        "results/deeptools/heatmaps/{heatmap}.mat.gz"
    output:
        "results/deeptools/heatmaps/{heatmap}.heatmap.pdf"
    params:
        lambda w: config["deeptools_params"][config["deeptools_heatmaps"][w.heatmap]["params"]]["heatmap"]
    log:
        "logs/deeptools/heatmaps/{heatmap}.plotHeatmap.log"
    shell:
        "plotHeatmap {params} -m {input} -out {output} 2> {log}"

