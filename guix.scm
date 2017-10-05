;;; PiGx_bsseq - reports pipeline for reads from bisulfite experiments.
;;; Copyright © 2017 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx_bsseq.
;;;
;;; PiGx_bsseq is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PiGx_bsseq:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; To install the package from a release tarball do this:
;;;
;;;  $ guix package --with-source=pigx_bsseq-0.0.1.tar.gz -f guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.13.0-2179-gbf3fa9965.

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix git-download)
             (guix build-system ant)
             (guix build-system gnu)
             (guix build-system r)
             (gnu packages)
             (gnu packages autotools)
             (gnu packages statistics)
             (gnu packages bioinformatics)
             (gnu packages compression)
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

(define-public r-methylkit-devel
  (let ((commit "46c8556f34eea9f068f2225e36adf11ba7ea6d07")
        (revision "1"))
    (package (inherit r-methylkit)
      (name "r-methylkit-devel")
      (version (string-append "1.3.3-" revision "." (string-take commit 9)))
      (source (origin
                (method git-fetch)
                (uri (git-reference
                      (url "https://github.com/al2na/methylKit.git")
                      (commit commit)))
                (file-name (string-append name "-" version "-checkout"))
                (sha256
                 (base32
                  "061yx5lgm5c37v9asnvbl4wxay04791cbxs52ar16x0a0gd13p53")))))))

(define-public r-annotationfilter
  (package
    (name "r-annotationfilter")
    (version "1.0.0")
    (source (origin
              (method url-fetch)
              (uri (bioconductor-uri "AnnotationFilter" version))
              (sha256
               (base32
                "0pxvswjzwibdfmrkdragxmzcl844z73pmkn82z92wahwa6gjfyi7"))))
    (properties
     `((upstream-name . "AnnotationFilter")))
    (build-system r-build-system)
    (propagated-inputs
     `(("r-genomicranges" ,r-genomicranges)
       ("r-lazyeval" ,r-lazyeval)))
    (home-page "https://github.com/Bioconductor/AnnotationFilter")
    (synopsis "Facilities for filtering Bioconductor annotation resources")
    (description
     "This package provides classes and other infrastructure to
implement filters for manipulating Bioconductor annotation resources.
The filters are used by ensembldb, Organism.dplyr, and other
packages.")
    (license artistic2.0)))

(define-public r-ensembldb
  (package
    (name "r-ensembldb")
    (version "2.0.4")
    (source
     (origin
       (method url-fetch)
       (uri (bioconductor-uri "ensembldb" version))
       (sha256
        (base32
         "1np96nry1hba8lk4bg3grf8w3k6xz9lgd2jcl3vrj6wsl184c3fr"))))
    (build-system r-build-system)
    (propagated-inputs
     `(("r-annotationdbi" ,r-annotationdbi)
       ("r-annotationfilter" ,r-annotationfilter)
       ("r-annotationhub" ,r-annotationhub)
       ("r-biobase" ,r-biobase)
       ("r-biocgenerics" ,r-biocgenerics)
       ("r-biostrings" ,r-biostrings)
       ("r-curl" ,r-curl)
       ("r-dbi" ,r-dbi)
       ("r-genomeinfodb" ,r-genomeinfodb)
       ("r-genomicfeatures" ,r-genomicfeatures)
       ("r-genomicranges" ,r-genomicranges)
       ("r-iranges" ,r-iranges)
       ("r-protgenerics" ,r-protgenerics)
       ("r-rsamtools" ,r-rsamtools)
       ("r-rsqlite" ,r-rsqlite)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-s4vectors" ,r-s4vectors)))
    (home-page "https://github.com/jotsetung/ensembldb")
    (synopsis "Utilities to create and use Ensembl-based annotation databases")
    (description
     "The package provides functions to create and use transcript
centric annotation databases/packages.  The annotation for the
databases are directly fetched from Ensembl using their Perl API.  The
functionality and data is similar to that of the TxDb packages from
the GenomicFeatures package, but, in addition to retrieve all
gene/transcript models and annotations from the database, the
ensembldb package provides also a filter framework allowing to
retrieve annotations for specific entries like genes encoded on a
chromosome region or transcript models of lincRNA genes.")
    ;; No version specified
    (license lgpl3+)))

(define-public r-organismdbi
  (package
    (name "r-organismdbi")
    (version "1.18.0")
    (source
     (origin
       (method url-fetch)
       (uri (bioconductor-uri "OrganismDbi" version))
       (sha256
        (base32
         "17jamgx9hqyi8ia48whqf7jj6ibdah2641zvx1xpv2lm8mhl3qzc"))))
    (properties `((upstream-name . "OrganismDbi")))
    (build-system r-build-system)
    (propagated-inputs
     `(("r-annotationdbi" ,r-annotationdbi)
       ("r-biobase" ,r-biobase)
       ("r-biocgenerics" ,r-biocgenerics)
       ("r-biocinstaller" ,r-biocinstaller)
       ("r-genomicfeatures" ,r-genomicfeatures)
       ("r-genomicranges" ,r-genomicranges)
       ("r-graph" ,r-graph)
       ("r-iranges" ,r-iranges)
       ("r-rbgl" ,r-rbgl)
       ("r-rsqlite" ,r-rsqlite)
       ("r-s4vectors" ,r-s4vectors)))
    (home-page "http://bioconductor.org/packages/OrganismDbi")
    (synopsis "Software to enable the smooth interfacing of database packages")
    (description "The package enables a simple unified interface to
several annotation packages each of which has its own schema by taking
advantage of the fact that each of these packages implements a select
methods.")
    (license artistic2.0)))

(define-public r-biovizbase
  (package
    (name "r-biovizbase")
    (version "1.24.0")
    (source
     (origin
       (method url-fetch)
       (uri (bioconductor-uri "biovizBase" version))
       (sha256
        (base32
         "1pfyhjwlxw9p2q5ip0irxpwndgakvn6z6ay5ahgz2gkkk8x8i29w"))))
    (properties `((upstream-name . "biovizBase")))
    (build-system r-build-system)
    (propagated-inputs
     `(("r-annotationdbi" ,r-annotationdbi)
       ("r-annotationfilter" ,r-annotationfilter)
       ("r-biocgenerics" ,r-biocgenerics)
       ("r-biostrings" ,r-biostrings)
       ("r-dichromat" ,r-dichromat)
       ("r-ensembldb" ,r-ensembldb)
       ("r-genomeinfodb" ,r-genomeinfodb)
       ("r-genomicalignments" ,r-genomicalignments)
       ("r-genomicfeatures" ,r-genomicfeatures)
       ("r-genomicranges" ,r-genomicranges)
       ("r-hmisc" ,r-hmisc)
       ("r-iranges" ,r-iranges)
       ("r-rcolorbrewer" ,r-rcolorbrewer)
       ("r-rsamtools" ,r-rsamtools)
       ("r-s4vectors" ,r-s4vectors)
       ("r-scales" ,r-scales)
       ("r-summarizedexperiment" ,r-summarizedexperiment)
       ("r-variantannotation" ,r-variantannotation)))
    (home-page "http://bioconductor.org/packages/biovizBase")
    (synopsis "Basic graphic utilities for visualization of genomic data")
    (description
     "The biovizBase package is designed to provide a set of
utilities, color schemes and conventions for genomic data.  It serves
as the base for various high-level packages for biological data
visualization.  This saves development effort and encourages
consistency.")
    (license artistic2.0)))

(define-public r-ggbio
  (package
    (name "r-ggbio")
    (version "1.24.1")
    (source
     (origin
       (method url-fetch)
       (uri (bioconductor-uri "ggbio" version))
       (sha256
        (base32
         "1nfqab0bvcs7rsh06qjj3i11hvz26gdrdc4yqbvs8b1x113kc0bs"))))
    (build-system r-build-system)
    (propagated-inputs
     `(("r-annotationdbi" ,r-annotationdbi)
       ("r-annotationfilter" ,r-annotationfilter)
       ("r-biobase" ,r-biobase)
       ("r-biocgenerics" ,r-biocgenerics)
       ("r-biostrings" ,r-biostrings)
       ("r-biovizbase" ,r-biovizbase)
       ("r-bsgenome" ,r-bsgenome)
       ("r-ensembldb" ,r-ensembldb)
       ("r-genomeinfodb" ,r-genomeinfodb)
       ("r-genomicalignments" ,r-genomicalignments)
       ("r-genomicfeatures" ,r-genomicfeatures)
       ("r-genomicranges" ,r-genomicranges)
       ("r-ggally" ,r-ggally)
       ("r-ggplot2" ,r-ggplot2)
       ("r-gridextra" ,r-gridextra)
       ("r-gtable" ,r-gtable)
       ("r-hmisc" ,r-hmisc)
       ("r-iranges" ,r-iranges)
       ("r-organismdbi" ,r-organismdbi)
       ("r-reshape2" ,r-reshape2)
       ("r-rsamtools" ,r-rsamtools)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-s4vectors" ,r-s4vectors)
       ("r-scales" ,r-scales)
       ("r-summarizedexperiment" ,r-summarizedexperiment)
       ("r-variantannotation" ,r-variantannotation)))
    (home-page "http://tengfei.github.com/ggbio/")
    (synopsis "Visualization tools for genomic data")
    (description
     "The ggbio package extends and specializes the grammar of
graphics for biological data.  The graphics are designed to answer
common scientific questions, in particular those often asked of high
throughput genomics data.  All core Bioconductor data structures are
supported, where appropriate.  The package supports detailed views of
particular genomic regions, as well as genome-wide overviews.
Supported overviews include ideograms and grand linear views.
High-level plots include sequence fragment length, edge-linked
interval to data view, mismatch pileup, and several splicing
summaries.")
    (license artistic2.0)))

(define %pigx-bsseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-bsseq
  (package
    (name "pigx_bsseq")
    (version %pigx-bsseq-version)
    (source (string-append (getcwd) "/pigx_bsseq-" version ".tar.gz"))
    (build-system gnu-build-system)
    (arguments
     `(#:tests? #f ; requires network access
       #:phases
       (modify-phases %standard-phases
         (add-after 'install 'wrap-executable
           ;; Make sure the executable finds all R modules.
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (wrap-program (string-append out "/bin/pigx-bsseq")
                 `("R_LIBS_SITE" ":" = (,(getenv "R_LIBS_SITE")))))
             #t)))))
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)))
    (inputs
     `(("r-minimal" ,r-minimal)
       ("r-annotationhub" ,r-annotationhub)
       ("r-dt" ,r-dt)
       ("r-genomation" ,r-genomation)
       ("r-genomeinfodb" ,r-genomeinfodb)
       ("r-methylkit-devel" ,r-methylkit-devel)
       ("r-rtracklayer" ,r-rtracklayer)
       ("r-rmarkdown" ,r-rmarkdown)
       ("r-bookdown" ,r-bookdown)
       ("r-ggplot2" ,r-ggplot2)
       ("r-ggbio" ,r-ggbio)
       ("ghc-pandoc" ,ghc-pandoc)
       ("ghc-pandoc-citeproc" ,ghc-pandoc-citeproc)
       ("python-wrapper" ,python-wrapper)
       ("python-pyyaml" ,python-pyyaml)
       ("snakemake" ,snakemake)
       ("bismark" ,bismark)
       ("fastqc" ,fastqc)
       ("bowtie" ,bowtie)
       ("trim-galore" ,trim-galore)
       ("cutadapt" ,cutadapt)
       ("samtools" ,samtools)))
    (home-page "https://github.com/BIMSBbioinfo/pigx_bsseq/")
    (synopsis "Bisulfite sequencing pipeline from fastq to methylation reports")
    (description "PiGx is a data processing pipeline for raw fastq
read data of bisulfite experiments; it produces reports on aggregate
methylation and coverage and can be used to produce information on
differential methylation and segmentation.")
    (license gpl3+)))

pigx-bsseq
