;;; PiGx_chipseq - reports pipeline for reads from ChIPseq experiments.
;;; Copyright © 2017 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx_chipseq.
;;;
;;; PiGx_chipseq is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PiGx_chipseq:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; To install the package from a release tarball do this:
;;;
;;;  $ guix package --with-source=pigx_chipseq-0.0.1.tar.gz -f guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.14.0-3177-gbcddf30af

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages autotools)
             (gnu packages statistics)
             (gnu packages bioinformatics)
             (gnu packages compression)
             (gnu packages check)
             (gnu packages cran)
             (gnu packages haskell)
             (gnu packages java)
             (gnu packages perl)
             (gnu packages python)
             (gnu packages web))

(define %pigx-chipseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-chipseq
  (package
    (name "pigx_chipseq")
    (version %pigx-chipseq-version)
    (source (string-append (getcwd) "/pigx_chipseq-" version ".tar.gz"))
    (build-system gnu-build-system)
    (arguments
     `(#:tests? #f ; requires network access
       #:phases
       (modify-phases %standard-phases
         (add-after 'install 'wrap-executable
           ;; Make sure the executable finds all R modules.
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (wrap-program (string-append out "/bin/pigx-chipseq")
                 `("R_LIBS_SITE" ":" = (,(getenv "R_LIBS_SITE")))
                 `("PYTHONPATH"  ":" = (,(getenv "PYTHONPATH")))))
             #t)))))
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)))
    (inputs
     `(("r-minimal" ,r-minimal)
       ("r-argparser" ,r-argparser)
       ("r-chipseq" ,r-chipseq)
       ("r-data-table" ,r-data-table)
       ("r-genomation" ,r-genomation)
       ("r-genomicranges" ,r-genomicranges)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-stringr" ,r-stringr)
       ("r-jsonlite" ,r-jsonlite)
       ("r-heatmaply" ,r-heatmaply)
       ("r-ggplot2" ,r-ggplot2)
       ("r-plotly" ,r-plotly)
       ("r-rmarkdown" ,r-rmarkdown)
       ("python-wrapper" ,python-wrapper)
       ("python-pyyaml" ,python-pyyaml)
       ("python-pytest" ,python-pytest)
       ("python-xlrd" ,python-xlrd)
       ("python-magic" ,python-magic)
       ("python-sh" ,python-sh)
       ("trim-galore" ,trim-galore)
       ("macs" ,macs)
       ("multiqc" ,multiqc)
       ("perl" ,perl)
       ("ghc-pandoc" ,ghc-pandoc)
       ("ghc-pandoc-citeproc" ,ghc-pandoc-citeproc)
       ("fastqc" ,fastqc)
       ("bowtie" ,bowtie)
       ("idr" ,idr)
       ("snakemake" ,snakemake)
       ("samtools" ,samtools)
       ("bedtools" ,bedtools)
       ("kentutils" ,kentutils)))
    (home-page "https://github.com/BIMSBbioinfo/pigx_chipseq/")
    (synopsis "TODO")
    (description "TODO")
    (license gpl3+)))

pigx-chipseq
