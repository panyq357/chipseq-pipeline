

ruleorder: bamCompare > bamCoverage


rule bamCoverage:
    input:
        bam = lambda w: config["bws"][w.sample]["b1"]
    output:
        bw = "results/BigWig/{sample}.bw"
    params:
        mappable_genome_size = config["mappable_genome_size"]
    threads:
        2
    resources:
        io = 100
    log:
        "logs/BigWig/{sample}.log"
    shell:
        '''
        bamCoverage \
            --numberOfProcessors {threads} \
            --bam {input.bam} \
            --outFileName {output.bw} \
            --effectiveGenomeSize {params.mappable_genome_size} \
            --binSize 20 \
            --normalizeUsing BPM \
            --smoothLength 60 \
            --extendReads \
            --centerReads \
            >> {log} 2>&1
        '''


rule bamCompare:
    input:
        b1 = lambda w: config["bws"][w.sample]["b1"],
        b2 = lambda w: config["bws"][w.sample]["b2"]
    output:
        bw = "results/BigWig/{sample}.bw"
    params:
        mappable_genome_size = config["mappable_genome_size"]
    threads:
        2
    resources:
        io = 100
    log:
        "logs/BigWig/{sample}.log"
    shell:
        '''
        bamCompare \
            --numberOfProcessors {threads} \
            --bamfile1 {input.b1} \
            --bamfile2 {input.b2} \
            --outFileName {output.bw} \
            --scaleFactorsMethod None \
            --binSize 20 \
            --normalizeUsing BPM \
            --smoothLength 60 \
            --extendReads \
            --centerReads \
            --outFileName {output} >> {log} 2>&1
        '''


rule computeMatrix:
    input:
        bw = lambda w: config["deeptools_heatmaps"][w.heatmap]["bw"],
        bed = lambda w: config["deeptools_heatmaps"][w.heatmap]["bed"]
    output:
        temp("results/deeptools_heatmaps/{heatmap}.mat.gz")
    params:
        lambda w: config["deeptools_params"][config["deeptools_heatmaps"][w.heatmap]["params"]]["matrix"]
    threads:
        20
    log:
        "logs/deeptools_heatmaps/{heatmap}.computeMatrix.log"
    shell:
        '''
        computeMatrix {params} \
            --numberOfProcessors {threads} \
            --regionsFileName {input.bed} \
            --scoreFileName {input.bw} \
            --outFileName {output}
        '''


rule plotHeatmap:
    input:
        "results/deeptools_heatmaps/{heatmap}.mat.gz"
    output:
        "results/deeptools_heatmaps/{heatmap}.heatmap.pdf"
    params:
        lambda w: config["deeptools_params"][config["deeptools_heatmaps"][w.heatmap]["params"]]["heatmap"]
    shell:
        '''
        plotHeatmap {params} -m {input} -out {output}
        '''

