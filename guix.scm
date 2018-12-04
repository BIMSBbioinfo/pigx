;;; PiGx_rnaseq - reports pipeline for RNAseq experiments.
;;; Copyright © 2017, 2018 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx_rnaseq.
;;;
;;; PiGx_rnaseq is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PiGx_rnaseq:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; To install the package from a release tarball do this:
;;;
;;;  $ guix package --with-source=pigx_rnaseq-0.0.1.tar.gz -f guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.15.0-4492-g91a4863d9

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix git-download)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages autotools)
             (gnu packages statistics)
             (gnu packages bioinformatics)
             (gnu packages compression)
             (gnu packages cran)
             (gnu packages haskell)
             (gnu packages java)
             (gnu packages python)
             (gnu packages web))

(define %pigx-rnaseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-rnaseq
  (package
    (name "pigx_rnaseq")
    (version %pigx-rnaseq-version)
    (source (string-append (getcwd) "/pigx_rnaseq-" version ".tar.gz"))
    (build-system gnu-build-system)
    (arguments '(#:parallel-tests? #f))  ; not supported
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)))
    (inputs
     `(("gzip" ,gzip)
       ("snakemake" ,snakemake)
       ("fastqc" ,fastqc)
       ("multiqc" ,multiqc)
       ("star" ,star)
       ("trim-galore" ,trim-galore)
       ("htseq" ,htseq)
       ("samtools" ,samtools)
       ("bedtools" ,bedtools)
       ("r-minimal" ,r-minimal)
       ("r-rmarkdown" ,r-rmarkdown)
       ("r-ggplot2" ,r-ggplot2)
       ("r-ggrepel" ,r-ggrepel)
       ("r-gprofiler" ,r-gprofiler)
       ("r-deseq2" ,r-deseq2)
       ("r-dt" ,r-dt)
       ("r-knitr" ,r-knitr)
       ("r-pheatmap" ,r-pheatmap)
       ("r-corrplot" ,r-corrplot)
       ("r-reshape2" ,r-reshape2)
       ("r-plotly" ,r-plotly)
       ("r-scales" ,r-scales)
       ("r-summarizedexperiment" ,r-summarizedexperiment)
       ("r-crosstalk" ,r-crosstalk)
       ("r-tximport" ,r-tximport)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-rjson" ,r-rjson)
       ("salmon" ,salmon)
       ("ghc-pandoc" ,ghc-pandoc)
       ("ghc-pandoc-citeproc" ,ghc-pandoc-citeproc)
       ("python-wrapper" ,python-wrapper)
       ("python-pyyaml" ,python-pyyaml)))
    (home-page "https://github.com/BIMSBbioinfo/pigx_rnaseq/")
    (synopsis "Analysis pipeline for single-cell RNA sequencing experiments")
    (description "PiGX scRNAseq is an analysis pipeline for preprocessing and
quality control for single cell RNA sequencing experiments.  The inputs are
read files from the sequencing experiment, and a configuration file which
describes the experiment.  It produces processed files for downstream analysis
and interactive quality reports.  The pipeline is designed to work with UMI
based methods.")
    (license gpl3+)))

pigx-rnaseq
