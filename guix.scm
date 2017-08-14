;;; PIGx_bsseq - reports pipeline for reads from bisulfite experiments.
;;; Copyright © 2017 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PIGx_bsseq.
;;;
;;; PIGx_bsseq is free software; see COPYING file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PIGx_bsseq:
;;;
;;;  $ export GUIX_PACKAGE_PATH=/path/to/guix-bimsb
;;;  $ guix environment -l guix.scm

(use-modules (guix packages)
             (guix licenses)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages statistics)
             (gnu packages bioinformatics)
             (gnu packages haskell)
             (gnu packages python)
             ;; These are modules from guix-bimsb
             ;; https://github.com/BIMSBbioinfo/guix-bimsb
             (bimsb packages staging)
             (bimsb packages bioinformatics-variants))

(package
  (name "pigx-bsseq")
  (version "0.0.0")
  (source #f)
  (build-system gnu-build-system)
  (inputs
   `(("r-minimal" ,r-minimal)
     ("r-annotationhub" ,r-annotationhub)
     ("r-dt" ,r-dt)
     ("r-genomation" ,r-genomation)
     ("r-methylkit-devel" ,r-methylkit-devel)
     ("r-rtracklayer" ,r-rtracklayer)
     ("r-rmarkdown" ,r-rmarkdown)
     ("r-bookdown" ,r-bookdown)
     ("ghc-pandoc" ,ghc-pandoc)
     ("ghc-pandoc-citeproc" ,ghc-pandoc-citeproc)
     ("python-wrapper" ,python-wrapper)
     ("snakemake" ,snakemake)
     ("bismark" ,bismark)
     ("fastqc" ,fastqc)
     ("bowtie" ,bowtie)
     ("trim-galore" ,trim-galore)
     ("cutadapt" ,cutadapt)
     ("samtools" ,samtools)))
  (home-page "TODO")
  (synopsis "TODO")
  (description "TODO")
  (license gpl3+))
