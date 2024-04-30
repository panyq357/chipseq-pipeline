configfile: "config/callpeak-config.yaml"

ruleorder: macs_with_input > macs_no_input


def get_bam(sample):
    return f"results/bwa_mapping/{sample}.bam"


rule macs_with_input:
    input:
        ip_bam = lambda w: get_bam(config["peaks"][w.peak]["ip"]),
        input_bam = lambda w: get_bam(config["peaks"][w.peak]["input"])
    output:
        peak = "results/macs_callpeak/{peak}/{peak}_peaks.xls"
    params:
        bam_format = config["macs_format"],
        mappable_genome_size = config["mappable_genome_size"],
        name = "{peak}",
        peak_type = lambda w: config["peaks"][w.peak]["type"]
    resources:
        disk_io = 5
    log:
        "logs/macs_callpeak/{peak}.log"
    run:
        if params.peak_type == "narrow":
            shell('''macs2 callpeak \
                        -t {input.ip_bam} \
                        -c {input.input_bam} \
                        -f {params.bam_format} \
                        -g {params.mappable_genome_size} \
                        -n {params.name} \
                        --outdir {output.out_dir} >> {log} 2>&1''')
        elif params.peak_type == "broad":
            shell('''macs2 callpeak \
                        -t {input.ip_bam} \
                        -c {input.input_bam} \
                        -f {params.bam_format} \
                        -g {params.mappable_genome_size} \
                        -n {params.name} \
                        --broad \
                        --outdir {output.out_dir} >> {log} 2>&1''')


rule macs_no_input:
    input:
        ip_bam = lambda w: get_bam(config["peaks"][w.peak]["ip"])
    output:
        peak = "results/macs_callpeak/{peak}/{peak}_peaks.xls"
    params:
        bam_format = config["macs_format"],
        mappable_genome_size = config["mappable_genome_size"],
        name = "{peak}",
        peak_type = lambda w: config["peaks"][w.peak]["type"]
    resources:
        disk_io = 5
    log:
        "logs/macs_call_peak/{peak}.log"
    run:
        if params.peak_type == "narrow":
            shell('''macs2 callpeak \
                        -t {input.ip_bam} \
                        -f {params.bam_format} \
                        -g {params.mappable_genome_size} \
                        -n {params.name} \
                        --outdir {output.out_dir} >> {log} 2>&1''')
        elif params.peak_type == "broad":
            shell('''macs2 callpeak \
                        -t {input.ip_bam} \
                        -f {params.bam_format} \
                        -g {params.mappable_genome_size} \
                        -n {params.name} \
                        --broad \
                        --outdir {output.out_dir} >> {log} 2>&1''')

