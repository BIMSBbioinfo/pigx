[![install with Guix badge](https://img.shields.io/badge/install%20with-guix-F4BB15.svg)](https://gnu.org/s/guix)

<div align="center">
<img src="logo.svg" alt="PiGx logo" width="30%" height="30%" />
</div>

# What is PiGx?

PiGx is a collection of genomics pipelines. More information can be found in [PiGx website](http://bioinformatics.mdc-berlin.de/pigx/)

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
pigx [pipeline] --init=settings
```

To generate a sample sheet template for any pipeline:

```sh
pigx [pipeline] --init=sample-sheet
```

Here's a simple example to run the RNAseq pipeline:

```sh
pigx rnaseq my-sample-sheet.csv --settings my-settings.yaml
```

To see all available options run `pigx --help`.


# Install

Pre-built binaries for PiGx are available through GNU Guix, the
functional package manager for reproducible, user-controlled software
management.  Install the complete pipeline bundle with the following
command:

```sh
guix package -i pigx
```

If you want to install PiGx from source, please make sure that all
required dependencies are installed and then follow the common GNU
build system steps after unpacking the [latest release
tarball](https://github.com/BIMSBbioinfo/pigx/releases/latest):

```sh
./configure --prefix=/some/where
make install
```

You can enable or disable each of the pipelines with the
`--enable-PIPELINE` and `--disable-PIPELINE` arguments to the
`configure` script.  `PIPELINE` is one of `bsseq`, `rnaseq`,
`scrnaseq`, `chipseq`, and `crispr`.  For more options run
`./configure --help`.


# License

All PiGx pipelines are free software: you can redistribute PiGx and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

See `LICENSE` for the full license text.
