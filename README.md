This is a pipeline for rice and arabidopsis ChIP-seq analysis.

## How to run this pipeline

Step 1. Filling `config/reads2bwa-config.yaml` with bwa index and raw FASTQ paths.

Step 2. Filling `config/callpeak-config.yaml` with peak's ip and input sample names
(must be matching sample names in `config/reads2bwa-config.yaml`).

Step 3. Filling `config/deeptools-config.yaml` with peak visualization config.

Step 4. Filling `config/binding-motif-config.yaml` for meme motif anlysis job info.

Step 5. Filling `config/diffbind-config.yaml` with DiffBind sample sheets,
for differential binding affinity analysis.

Step 6. Filling `config/peak-to-gene-config.yaml` with GTF path for assign peaks
to its overlapping gene.

Step 7. For arabidopsis annotation, uncomment these lines in `workflow/Snakefile`.

```
include: "rules/peak-anno-at.smk"
```

```
# peak-anno-at
expand("results/peak_anno/{job}/{job}.at_anno_diffbind.xlsx", job=config["diffbind-sample-sheets"]),
expand("results/peak_anno/{peak}/{peak}.at_anno_macs.xlsx", peak=config["peaks"])
```

Finally, run this.

```bash
snakemake --cores 20 --resources disk_io=10
```

## Results explained

```txt
results/
├── bamCompare                                                           # bigwig files for IP minus Input.
│   ├── JA-1-compare.bw
│   ├── JA-2-compare.bw
│   ├── JA-3-compare.bw
│   ├── Mock-1-compare.bw
│   ├── Mock-2-compare.bw
│   └── Mock-3-compare.bw
├── bamCoverage                                                          # coverage of BAms.
│   ├── Input_JA-1.bw
│   ├── Input_JA-2.bw
│   ├── Input_JA-3.bw
│   ├── Input_Mock-1.bw
│   ├── Input_Mock-2.bw
│   ├── Input_Mock-3.bw
│   ├── IP_JA-1.bw
│   ├── IP_JA-2.bw
│   ├── IP_JA-3.bw
│   ├── IP_Mock-1.bw
│   ├── IP_Mock-2.bw
│   └── IP_Mock-3.bw
├── bwa_mapping                                                          # BWA mapping results.
│   ├── Input_JA-1.bam
│   ├── Input_JA-1.bam.bai
│   ├── Input_JA-1.flagstats.txt
│   ├── Input_JA-2.bam
│   ├── Input_JA-2.bam.bai
│   ├── Input_JA-2.flagstats.txt
│   ├── Input_JA-3.bam
│   ├── Input_JA-3.bam.bai
│   ├── Input_JA-3.flagstats.txt
│   ├── Input_Mock-1.bam
│   ├── Input_Mock-1.bam.bai
│   ├── Input_Mock-1.flagstats.txt
│   ├── Input_Mock-2.bam
│   ├── Input_Mock-2.bam.bai
│   ├── Input_Mock-2.flagstats.txt
│   ├── Input_Mock-3.bam
│   ├── Input_Mock-3.bam.bai
│   ├── Input_Mock-3.flagstats.txt
│   ├── IP_JA-1.bam
│   ├── IP_JA-1.bam.bai
│   ├── IP_JA-1.flagstats.txt
│   ├── IP_JA-2.bam
│   ├── IP_JA-2.bam.bai
│   ├── IP_JA-2.flagstats.txt
│   ├── IP_JA-3.bam
│   ├── IP_JA-3.bam.bai
│   ├── IP_JA-3.flagstats.txt
│   ├── IP_Mock-1.bam
│   ├── IP_Mock-1.bam.bai
│   ├── IP_Mock-1.flagstats.txt
│   ├── IP_Mock-2.bam
│   ├── IP_Mock-2.bam.bai
│   ├── IP_Mock-2.flagstats.txt
│   ├── IP_Mock-3.bam
│   ├── IP_Mock-3.bam.bai
│   └── IP_Mock-3.flagstats.txt
├── diffbind_analyze                                                     # DiffBind analysis results.
│   ├── no-input
│   └── with-input
├── diffbind_counts                                                      # DiffBind counts matrix.
│   ├── no-input_counts.rds
│   ├── no-input_peakset.csv
│   ├── with-input_counts.rds
│   └── with-input_peakset.csv
├── diff_peak_heatmap                                                    # Heatmap of Differential binding peaks detected by DiffBind.
│   ├── no-input.diff_peak_heatmap.png
│   └── with-input.diff_peak_heatmap.png
├── fastp                                                                # raw FASTQ QC reports.
│   ├── Input_JA-1.html
│   ├── Input_JA-1.json
│   ├── Input_JA-2.html
│   ├── Input_JA-2.json
│   ├── Input_JA-3.html
│   ├── Input_JA-3.json
│   ├── Input_Mock-1.html
│   ├── Input_Mock-1.json
│   ├── Input_Mock-2.html
│   ├── Input_Mock-2.json
│   ├── Input_Mock-3.html
│   ├── Input_Mock-3.json
│   ├── IP_JA-1.html
│   ├── IP_JA-1.json
│   ├── IP_JA-2.html
│   ├── IP_JA-2.json
│   ├── IP_JA-3.html
│   ├── IP_JA-3.json
│   ├── IP_Mock-1.html
│   ├── IP_Mock-1.json
│   ├── IP_Mock-2.html
│   ├── IP_Mock-2.json
│   ├── IP_Mock-3.html
│   └── IP_Mock-3.json
├── get_peak_fasta                                                       # FASTA for meme-chip motif analysis.
│   ├── Mock-1-with-input-top-1000-radius-50.fasta
│   ├── Mock-2-with-input-top-1000-radius-50.fasta
│   └── Mock-3-with-input-top-1000-radius-50.fasta
├── macs_callpeak                                                        # MACS callpeak results.
│   ├── JA-1-no-input
│   ├── JA-1-with-input
│   ├── JA-2-no-input
│   ├── JA-2-with-input
│   ├── JA-3-no-input
│   ├── JA-3-with-input
│   ├── Mock-1-no-input
│   ├── Mock-1-with-input
│   ├── Mock-2-no-input
│   ├── Mock-2-with-input
│   ├── Mock-3-no-input
│   └── Mock-3-with-input
├── meme                                                                 # meme-chip results.
│   ├── Mock-1-with-input-top-1000-radius-50
│   ├── Mock-2-with-input-top-1000-radius-50
│   └── Mock-3-with-input-top-1000-radius-50
├── peak_anno                                                            # Peak annotation results.
│   ├── JA-1-no-input
│   ├── JA-1-with-input
│   ├── JA-2-no-input
│   ├── JA-2-with-input
│   ├── JA-3-no-input
│   ├── JA-3-with-input
│   ├── Mock-1-no-input
│   ├── Mock-1-with-input
│   ├── Mock-2-no-input
│   ├── Mock-2-with-input
│   ├── Mock-3-no-input
│   ├── Mock-3-with-input
│   ├── no-input
│   └── with-input
├── peak_heatmap                                                         # Peak heatmap for QC.
│   ├── JA-1-no-input.heatmap.png
│   ├── JA-1-with-input.heatmap.png
│   ├── JA-2-no-input.heatmap.png
│   ├── JA-2-with-input.heatmap.png
│   ├── JA-3-no-input.heatmap.png
│   ├── JA-3-with-input.heatmap.png
│   ├── Mock-1-no-input.heatmap.png
│   ├── Mock-1-with-input.heatmap.png
│   ├── Mock-2-no-input.heatmap.png
│   ├── Mock-2-with-input.heatmap.png
│   ├── Mock-3-no-input.heatmap.png
│   └── Mock-3-with-input.heatmap.png
└── tss_to_tes_heatmap                                                   # TSS to TES heatmap for viewing peak location.
    ├── compare.tss_to_tes.heatmap.png
    ├── compare.tss_to_tes.profile.png
    ├── coverage.tss_to_tes.heatmap.png
    └── coverage.tss_to_tes.profile.png
```

