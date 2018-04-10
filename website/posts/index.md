title: PiGx: Pipelines in Genomics
---

# What is PiGx?

PiGx is a collection of genomics pipelines.  All pipelines are easily
configured with a simple sample sheet and a descriptive settings file.
The result is a set of comprehensive, interactive HTML reports with
interesting findings about your samples.

###### BUTTONS HERE

PiGx includes the following pipelines:

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

Here are some teaser snapshots taken from the HTML reports:

###### INSERT EXAMPLES HERE

# Getting started

To run PiGx on your experimental data you only need to follow these
three steps:

- **Describe your samples.** Describe all your samples in a CSV file
   `sample_sheet.csv`.  To generate a sample sheet template for any
   pipeline run this command:

       pigx [pipeline] --init=sample-sheet

- **Tweak the default settings.** Specify input and output
   directories, and override defaults by providing a `settings.yaml`.
   To generate a template for this file for any pipeline run this
   command:

       pigx [pipeline] --init=settings

- **Run the pipeline!** Here's a simple example to run the RNAseq
   pipeline when the sample sheet `sample_sheet.csv` is in the current
   working directory:

       pigx rnaseq -s settings.yaml

That's it!  After some time you'll get a bunch of reports that you can
view in your browser.

For more detailed information about each of the pipelines [see the
online documentation](http://bioinformatics.mdc-berlin.de/pigx_docs).


# Get it

Pre-built binaries for PiGx are available through [GNU
Guix](https://gnu.org/software/guix), the functional package manager
for reproducible, user-controlled software management.  Install the
complete pipeline bundle with the following command:

```sh
guix package -i pigx
```

PiGx is [free
software](https://www.fsf.org/about/what-is-free-software) under the
[GNU General Public License](https://www.gnu.org/licenses/gpl.html)
(version 3 or, at your option, any later version).  You can [get the
complete source code here](https://github.com/BIMSBbioinfo/pigx).
Your contributions are welcome!

Consider [subscribing here to infrequent announcements regarding
PiGx](https://groups.google.com/forum/#!forum/pigx-announcements/join),
such as release announcements.
