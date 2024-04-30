This is a pipeline for rice and arabidopsis ChIP-seq analysis.

## How to run this pipeline

Step 1. Filling config files in `config` directory.

Step 2. Comment and uncomment lines in `workflow/Snakefile` to specify what results need to be generated.

Finally, run this.

```bash
snakemake --cores 20 --resources disk_io=10
```

