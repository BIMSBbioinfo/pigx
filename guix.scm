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
;;; v0.13.0-4757-ga0837294a

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix build-system ant)
             (guix build-system gnu)
             (guix build-system r)
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

;; FIXME: This package includes pre-built Java classes.
(define-public fastqc
  (package
    (name "fastqc")
    (version "0.11.5")
    (source
     (origin
       (method url-fetch)
       (uri (string-append "http://www.bioinformatics.babraham.ac.uk/"
                           "projects/fastqc/fastqc_v"
                           version "_source.zip"))
       (sha256
        (base32
         "18rrlkhcrxvvvlapch4dpj6xc6mpayzys8qfppybi8jrpgx5cc5f"))))
    (build-system ant-build-system)
    (arguments
     `(#:tests? #f ; there are no tests
       #:build-target "build"
       #:phases
       (modify-phases %standard-phases
         ;; There is no installation target
         (replace 'install
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let* ((out   (assoc-ref outputs "out"))
                    (bin   (string-append out "/bin"))
                    (share (string-append out "/share/fastqc/"))
                    (exe   (string-append share "/fastqc")))
               (for-each mkdir-p (list bin share))
               (copy-recursively "bin" share)
               (substitute* exe
                 (("my \\$java_bin = 'java';")
                  (string-append "my $java_bin = '"
                                 (assoc-ref inputs "java")
                                 "/bin/java';")))
               (chmod exe #o555)
               (symlink exe (string-append bin "/fastqc"))
               #t))))))
    (inputs
     `(("java" ,icedtea)
       ("perl" ,perl)))  ; needed for the wrapper script
    (native-inputs
     `(("unzip" ,unzip)))
    (home-page "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/")
    (synopsis "Quality control tool for high throughput sequence data")
    (description
     "FastQC aims to provide a simple way to do some quality control
checks on raw sequence data coming from high throughput sequencing
pipelines.  It provides a modular set of analyses which you can use to
give a quick impression of whether your data has any problems of which
you should be aware before doing any further analysis.

The main functions of FastQC are:

@itemize
@item Import of data from BAM, SAM or FastQ files (any variant);
@item Providing a quick overview to tell you in which areas there may
  be problems;
@item Summary graphs and tables to quickly assess your data;
@item Export of results to an HTML based permanent report;
@item Offline operation to allow automated generation of reports
  without running the interactive application.
@end itemize\n")
    (license gpl3+)))

(define-public r-argparser
  (package
    (name "r-argparser")
    (version "0.4")
    (source
     (origin
       (method url-fetch)
       (uri (cran-uri "argparser" version))
       (sha256
        (base32
         "0s1wxshx4jk69wfxhycx973q6y8cmqrfymyjklhq1i8xrj0kmmx9"))))
    (build-system r-build-system)
    (home-page "https://bitbucket.org/djhshih/argparser")
    (synopsis "Command-line argument parser")
    (description
     "This package provides a cross-platform command-line argument parser
written purely in R with no external dependencies.  It is useful with the
Rscript front-end and facilitates turning an R script into an executable
script.")
     gpl3+)))

;; FIXME: has to be removed after guixr pull
(define-public r-heatmaply
 (package
   (name "r-heatmaply")
   (version "0.14.1")
   (source
     (origin
       (method url-fetch)
       (uri (cran-uri "heatmaply" version))
       (sha256
         (base32
           "03p2caclhfgqgpx3wwck5h06jy3mxgs05gjmwkb7hmwghkjh41jc"))))
   (build-system r-build-system)
   (propagated-inputs
     `(("r-assertthat" ,r-assertthat)
       ("r-colorspace" ,r-colorspace)
       ("r-dendextend" ,r-dendextend)
       ("r-ggplot2" ,r-ggplot2)
       ("r-gplots" ,r-gplots)
       ("r-htmlwidgets" ,r-htmlwidgets)
       ("r-magrittr" ,r-magrittr)
       ("r-plotly" ,r-plotly)
       ("r-rcolorbrewer" ,r-rcolorbrewer)
       ("r-reshape2" ,r-reshape2)
       ("r-scales" ,r-scales)
       ("r-seriation" ,r-seriation)
       ("r-viridis" ,r-viridis)
       ("r-webshot" ,r-webshot)))
   (home-page
     "https://cran.r-project.org/package=heatmaply")
   (synopsis
     "Interactive Cluster Heat Maps Using 'plotly'")
   (description
     "Create interactive cluster 'heatmaps' that can be saved as a stand- alone HTML file, embedded in 'R Markdown' documents or in a 'Shiny' app, and available in the 'RStudio' viewer pane.  Hover the mouse pointer over a cell to show details or drag a rectangle to zoom.  A 'heatmap' is a popular graphical method for visualizing high-dimensional data, in which a table of numbers are encoded as a grid of colored cells.  The rows and columns of the matrix are ordered to highlight patterns and are often accompanied by 'dendrograms'. 'Heatmaps' are used in many fields for visualizing observations, correlations, missing values patterns, and more.  Interactive 'heatmaps' allow the inspection of specific value by hovering the mouse over a cell, as well as zooming into a region of the 'heatmap' by dragging a rectangle around the relevant area.  This work is based on the 'ggplot2' and 'plotly.js' engine.  It produces similar 'heatmaps' as 'heatmap.2' or 'd3heatmap', with the advantage of speed ('plotly.js' is able to handle larger size matrix), the ability to zoom from the 'dendrogram' panes, and the placing of factor variables in the sides of the 'heatmap'.")
   (license #f)))

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
       ("r-chipseq", r-chipseq)
       ("r-data-table" ,r-data-table)
       ("r-genomation" ,r-genomation)
       ("r-genomicranges" ,r-genomicranges)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-rcas" ,r-rcas)
       ("r-stringr" ,r-stringr)
       ("r-jsonlite" ,r-jsonlite)
       ("r-heatmaply" ,r-heatmaply)
       ("r-ggplot2" ,r-ggplot2)
       ("r-plotly" ,r-plotly)
       ("python-wrapper" ,python-wrapper)
       ("python-pyyaml" ,python-pyyaml)
       ("python-pytest" ,python-pytest)
       ("snakemake" ,snakemake)
       ("macs" ,macs)
       ("multiqc", multiqc)
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
