title: PiGx: Supplementary Materials
---

These are supplementary files for the PiGx paper.


## PiGx scRNA-seq

The HTML report output of the scRNA-seq pipeline from the analysis of
the single cell sequencing dataset from (Hu et al. 2017) [can be
accessed
here](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/scrnaseq/mm10.scRNA-Seq.report.html).


## PiGx ChIP-seq

The HTML report output of the ChIP-seq pipeline from the analysis of
the ChIP-seq datasets from Hon. et al [can be accessed here](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/chipseq/ChIP_Seq_Report.html).


## PiGx RNA-seq

The supplementary results for the RNA-seq pipeline from the analysis
of the RNA-seq datasets from Hon. et al consist of seven zipped
folders. Each folder contains the differential expression analysis
results as listed in Table S1. For each differential expression
analysis folder, there are three distinct HTML report files (files
ending with `.deseq.report.html`) generated using

- STAR-based gene counts,

- Salmon-based gene counts,

- Salmon-based transcript counts.

For each HTML file, there are two supporting tab-separated
files:

- The full table of DESeq2 differential expression statistics for each
  gene/transcript (files ending with `.deseq_results.tsv`),

- The full table of normalized counts (using DESeq2 variance
  stabilizing transformation) (files ending with
  `.normalized_counts.tsv`).

Below are links for the supplementary files for each differential
expression analysis performed for the use-case RNA-seq datasets:

- [tet2_diff_day3](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_diff_day3.tgz)

- [tet2_diff_day6](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_diff_day6.tgz)

- [tet2_vs_WT_day0](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_vs_WT_day0.tgz)

- [tet2_vs_WT_day3](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_vs_WT_day3.tgz)

- [tet2_vs_WT_day6](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/tet2_vs_WT_day6.tgz)

- [WT_diff_day3](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/WT_diff_day3.tgz)

- [WT_diff_day6](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/rnaseq/WT_diff_day6.tgz)


## PiGx BS-seq

The supplementary results for the BS-seq pipeline from the analysis of
the RNA-seq datasets from Hon. et al consist of 9 reports in HTML
format; We have highlighted the most important figures from the second
biological replicate in the study (as used in the original study), and
for completeness, we include here the same reports for the first
biological replicate, as well as the post-maturation data from
labelled “day 3”.

Each HTML document represents a stand-alone summary of some the key
findings for the given sample, along with a description of the key
operations carried out by the pipeline.

PiGx BS-seq Analysis Title
Download Link
[Tet2 deletion (biological replicate 1) report](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/tet2_se_bt2.sorted.deduped_mm10_canon_final.html)
[Wild-type (biological replicate 1) report](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/WT_se_bt2.sorted.deduped_mm10_canon_final.html)
[Differential methylation, biological replicate 1, (tet2-/-)  - WT](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/diffmeth-report.0vs1.html)
[Tet2 deletion (biological replicate 2) report](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/tet2-brep2_se_bt2.sorted.deduped_mm10_canon_final.html)
[Wild type (biological replicate 2) report](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/WT-brep2_se_bt2.sorted.deduped_mm10_canon_final.html)
[Differential methylation, biological replicate 2, (tet2-/-)  - WT](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/diffmeth-report.4vs5.html)
[Tet2 deletion (day 3) report](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/tet2-day3_se_bt2.sorted.deduped_mm10_canon_final.html)
[Wild type (day 3) report](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/WT-day3_se_bt2.sorted.deduped_mm10_canon_final.html)
[Differential methylation, day 3, (tet2-/-)  - WT](http://bimsbstatic.mdc-berlin.de/akalin/PiGx/supplementary_material/bsseq/diffmeth-report.2vs3.html)


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