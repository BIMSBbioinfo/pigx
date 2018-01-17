title: PiGx: Pipelines in Genomics
---

# What is PiGx?

PiGx is a collection of genomics pipelines.

It includes the following pipelines:

- [PiGx BSseq](https://github.com/BIMSBbioinfo/pigx_bsseq) for raw
  fastq read data of bisulfite experiments

- [PiGx RNAseq](https://github.com/BIMSBbioinfo/pigx_rnaseq) for RNAseq samples

- [PiGx scRNAseq](https://github.com/BIMSBbioinfo/pigx_scrnaseq) for
  single cell dropseq analysis

- [PiGx ChIPseq](https://github.com/BIMSBbioinfo/pigx_chipseq) for
  reads from ChIPseq experiments

- [PiGx CRISPR](https://github.com/BIMSBbioinfo/pigx_crispr) *(work in progress)*
  for the analysis of sequence mutations in CRISPR-CAS9 targeted
  amplicon sequencing data

All pipelines are easily configured with a sample sheet (in CSV
format) and a descriptive settings file (in YAML format).  For more
detailed information see the README.md file for each of the pipelines
in the `pipelines` directory.


# Getting started

To run PiGx on your experimental data, describe your samples in a CSV
file `sample_sheet.csv`, provide a `settings.yaml` to override the
defaults defaults, and select the pipeline.

To generate a settings file template for any pipeline:

```sh
pigx [pipeline] --init-settings my-settings.yaml
```

To generate a sample sheet template for any pipeline:

```sh
pigx [pipeline] --init-sample-sheet my-sample-sheet.csv
```

Here's a simple example to run the RNAseq pipeline:

```sh
pigx rnaseq my-sample-sheet.csv --settings my-settings.yaml
```

To see all available options run `pigx --help`.


# Install

Pre-built binaries for PiGx are (soon) available through GNU Guix, the
functional package manager for reproducible, user-controlled software
management.  Install the complete pipeline bundle with the following
command: 

```sh
guix package -i pigx
```
