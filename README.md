
This a [Snakemake](https://snakemake.github.io/) pipeline for performing ChIP-seq data analysis.

Analysis steps can be divided into two type: jobs and modifications.

For jobs, input files and parameters must be explicitly specified in the configuration file in order for them to run.

For modifications, they are run implicitly by adding suffix/infix in filenames.

To run this pipeline, first, copy config and Snakefile.

```bash
cp -r config_template/ config
cp workflow/Snakefile_template workflow/Snakefile
```

After filling config files, run this pipeline using following command.

```bash
snakemake --cores 20 --resources io=100
```
