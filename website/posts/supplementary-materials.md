title: Supplementary Materials
---

These are supplementary files for the PiGx paper.  All reports were
generated with PiGx 0.0.3 as installed with [GNU
Guix](https://gnu.org/software/guix) at version
v0.14.0-7054-g5149aeb7e (git commit
5149aeb7e62cf62398b55be38469cd28c25d8d7d).

We provide a [generated Docker image on
Dockerhub](https://hub.docker.com/r/bimsbbioinfo/pigx/).  The tag
`bimsbbioinfo/pigx:publication` corresponds to the image at the above
Guix version, which was used for the publication.

We also provide a [SquashFS file system
image](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/pigx-singularity-image.squashfs),
which can be used with systems like Singularity.

The SquashFS image was created with this command:

```
guix pack pigx            \
          -f squashfs     \
          -S /bin=bin     \
          -S /lib=lib     \
          -S /share=share \
          glibc-utf8-locales tzdata coreutils bash
```

This command was used for generating the Docker image:

```
guix pack pigx            \
          -f docker       \
          -C none         \
          -S /bin=bin     \
          -S /lib=lib     \
          -S /share=share \
          glibc-utf8-locales tzdata coreutils bash
```


## PiGx scRNA-seq

The HTML report output of the scRNA-seq pipeline from the analysis of
the single cell sequencing dataset from (Hu et al. 2017) [can be
accessed
here](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/scrnaseq/mm10.scRNA-Seq.report.html).

An archive containing a sample sheet, settings file, and instructions
on how to reproduce the report [can be downloaded
here](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/scrnaseq/scrnaseq-notes.tar.gz).

## PiGx ChIP-seq

The HTML report output of the ChIP-seq pipeline from the analysis of
the ChIP-seq datasets from Hon. et al [can be accessed here](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/chipseq/ChIP_Seq_Report.html).

An archive containing a sample sheet, settings file, and instructions
on how to reproduce the report [can be downloaded
here](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/chipseq/chipseq-notes.tar.gz).


## PiGx RNA-seq

Below are links for the HTML reports for each differential expression 
analysis performed using STAR-based gene counts for the use-case 
RNA-seq datasets:

- [tet2_diff_day3](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_diff_day3.html)

- [tet2_diff_day6](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_diff_day6.html)

- [tet2_vs_WT_day0](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_vs_WT_day0.html)

- [tet2_vs_WT_day3](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_vs_WT_day3.html)

- [tet2_vs_WT_day6](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_vs_WT_day6.html)

- [WT_diff_day3](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/WT_diff_day3.html)

- [WT_diff_day6](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/WT_diff_day6.html)

An archive containing a sample sheet, settings file, and instructions
on how to reproduce the report [can be downloaded
here](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/rnaseq-notes.tar.gz).


## PiGx BS-seq

The supplementary results for the BS-seq pipeline from the analysis of
the RNA-seq datasets from Hon. et al consist of three reports in HTML
format; We have highlighted the most important figures from the study.

Each HTML document represents a stand-alone summary of some the key
findings for the given sample, along with a description of the key
operations carried out by the pipeline.

- [Tet2 deletion (biological replicate 1) report](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/tet2_se_bt2.sorted.deduped_mm10_final.html)

- [Wild-type (biological replicate 1) report](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/WT_se_bt2.sorted.deduped_mm10_final.html)

- [Differential methylation, biological replicate 1, (tet2-/-) - WT](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/diffmeth-report.0vs1.html)

An archive containing a sample sheet, settings file, and instructions
on how to reproduce the report [can be downloaded
here](https://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/bsseq-notes.tar.gz).


# About PiGx

PiGx is [free
software](https://www.fsf.org/about/what-is-free-software) under the
[GNU General Public License](https://www.gnu.org/licenses/gpl.html)
(version 3 or, at your option, any later version).  You can [get the
complete source code here](https://github.com/BIMSBbioinfo/pigx).
Your contributions are welcome!

Consider [subscribing here to infrequent announcements regarding
PiGx](https://groups.google.com/forum/#!forum/pigx-announcements/join),
such as release announcements.
