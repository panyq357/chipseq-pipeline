A pipeline for ChIP-seq analysis.

## How to run this pipeline

First, run these.

```bash
cp -r config_template config 
cp workflow/Snakefile_template workflow/Snakefile
```

Then filling config files in `config` directory, and specify the outputs in `workflow/Snakefile`.

Finnaly, run this.

```bash
snakemake --cores 20 --resources io=100
```

