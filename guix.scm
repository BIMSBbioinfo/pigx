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
;;;  $ guix environment -l guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.13.0-2179-gbf3fa9965.

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix git-download)
             (guix build-system ant)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages statistics)
             (gnu packages bioinformatics)
             (gnu packages compression)
             (gnu packages haskell)
             (gnu packages java)
             (gnu packages perl)
             (gnu packages python))

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
     ("r-genomeinfodb" ,r-genomeinfodb)
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
