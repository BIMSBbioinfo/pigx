;;; PiGx_scrnaseq - reports pipeline for single cell RNAseq experiments.
;;; Copyright © 2018, 2019 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx_scrnaseq.
;;;
;;; PiGx_scrnaseq is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PiGx_scrnaseq:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; To install the package from a release tarball do this:
;;;
;;;  $ guix package --with-source=pigx_scrnaseq-0.0.1.tar.gz -f guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.14.0-3089-g84c195e50

(use-modules (guix packages) (gnu packages))

(define %pigx-scrnaseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-scrnaseq-development
  (package (inherit (specification->package "pigx-scrnaseq"))
    (version %pigx-scrnaseq-version)
    (source (string-append (getcwd) "/pigx_scrnaseq-" version ".tar.gz"))
    (native-inputs
     `(("autoconf" ,(specification->package "autoconf"))
       ("automake" ,(specification->package "automake"))))))
    (inputs
     `(("coreutils" ,coreutils)
       ("perl" ,perl)
       ("dropseq-tools" ,dropseq-tools)
       ("fastqc" ,fastqc)
       ("java-picard" ,java-picard)
       ("java" ,icedtea-8)
       ("samtools", samtools)
       ("python-wrapper" ,python-wrapper)
       ("python-pyyaml" ,python-pyyaml)
       ("python-pandas" ,python-pandas)
       ("python-numpy" ,python-numpy)
       ("python-loompy" ,python-loompy)
       ("python-magic", python-magic)
       ("ghc-pandoc" ,ghc-pandoc-1)
       ("ghc-pandoc-citeproc" ,ghc-pandoc-citeproc-with-pandoc-1)
       ("snakemake" ,snakemake)
       ("star" ,star-2.5.2)
       ("r-minimal" ,r-minimal)
       ("r-argparser" ,r-argparser)
       ("r-cowplot" ,r-cowplot)
       ("r-data-table" ,r-data-table)
       ("r-delayedarray" ,r-delayedarray)
       ("r-delayedmatrixstats" ,r-delayedmatrixstats)
       ("r-dplyr" ,r-dplyr)
       ("r-dropbead" ,r-dropbead)
       ("r-dt" ,r-dt)
       ("r-genomicalignments" ,r-genomicalignments)
       ("r-genomicfiles" ,r-genomicfiles)
       ("r-genomicranges" ,r-genomicranges)
       ("r-ggplot2" ,r-ggplot2)
       ("r-hdf5array" ,r-hdf5array)
       ("r-pheatmap" ,r-pheatmap)
       ("r-rmarkdown" ,r-rmarkdown)
       ("r-rsamtools" ,r-rsamtools)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-rtsne" ,r-rtsne)
       ("r-scater" ,r-scater)
       ("r-scran" ,r-scran)
       ("r-singlecellexperiment" ,r-singlecellexperiment)
       ("r-stringr" ,r-stringr)
       ("r-seurat", r-seurat)
       ("r-yaml" ,r-yaml)))
    (home-page "https://github.com/BIMSBbioinfo/pigx_scrnaseq/")
    (synopsis "Analysis pipeline for single-cell RNA sequencing experiments")
    (description "PiGX scRNAseq is an analysis pipeline for preprocessing and quality control for single cell
RNA sequencing experiments.  The inputs are read files from the sequencing experiment, and a configuration
file which describes the experiment.  It produces processed files for downstream analysis and interactive
quality reports.  The pipeline is designed to work with UMI based methods.")
    (license gpl3+)))
