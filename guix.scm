;;; PiGx_bsseq - reports pipeline for reads from bisulfite experiments.
;;; Copyright © 2017, 2018 Ricardo Wurmus <rekado@elephly.net>
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
;;; TODO

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages autotools)
             (gnu packages base)
             (gnu packages statistics)
             (gnu packages bioinformatics)
             (gnu packages compression)
             (gnu packages cran)
             (gnu packages curl)
             (gnu packages haskell)
             (gnu packages java)
             (gnu packages ncurses)
             (gnu packages perl)
             (gnu packages pkg-config)
             (gnu packages python)
             (gnu packages tls)
             (gnu packages web))

(define-public my/htslib
  (package
    (name "htslib")
    (version "1.7")
    (source (origin
              (method url-fetch)
              (uri (string-append
                    "https://github.com/samtools/htslib/releases/download/"
                    version "/htslib-" version ".tar.bz2"))
              (sha256
               (base32
                "1il6i2p84b0y9c93dhvzzki1ifw9bvapm2mvpr0xvb2nq8jlwgdy"))))
    (build-system gnu-build-system)
    (inputs
     `(("openssl" ,openssl)
       ("curl" ,curl)
       ("zlib" ,zlib)))
    (native-inputs
     `(("perl" ,perl)))
    (home-page "http://www.htslib.org")
    (synopsis "C library for reading/writing high-throughput sequencing data")
    (description
     "HTSlib is a C library for reading/writing high-throughput sequencing
data.  It also provides the @command{bgzip}, @command{htsfile}, and
@command{tabix} utilities.")
    ;; Files under cram/ are released under the modified BSD license;
    ;; the rest is released under the Expat license
    (license (list expat bsd-3))))

(define-public my/samtools
  (package
    (name "samtools")
    (version "1.7")
    (source
     (origin
       (method url-fetch)
       (uri
        (string-append "mirror://sourceforge/samtools/samtools/"
                       version "/samtools-" version ".tar.bz2"))
       (sha256
        (base32
         "18acyqysbxpydlc44lqv2hpp57l06bs9a3yqmcvjk8va2xrrdc77"))))
    (build-system gnu-build-system)
    (arguments
     `(#:modules ((ice-9 ftw)
                  (ice-9 regex)
                  (guix build gnu-build-system)
                  (guix build utils))
       #:make-flags (list (string-append "prefix=" (assoc-ref %outputs "out")))
       #:configure-flags (list "--with-ncurses" "--with-htslib=system")
       #:phases
       (modify-phases %standard-phases
         (add-after 'unpack 'patch-tests
           (lambda _
             (substitute* "test/test.pl"
               ;; The test script calls out to /bin/bash
               (("/bin/bash") (which "bash")))
             #t))
         (add-after 'install 'install-library
           (lambda* (#:key outputs #:allow-other-keys)
             (let ((lib (string-append (assoc-ref outputs "out") "/lib")))
               (install-file "libbam.a" lib)
               #t)))
         (add-after 'install 'install-headers
           (lambda* (#:key outputs #:allow-other-keys)
             (let ((include (string-append (assoc-ref outputs "out")
                                           "/include/samtools/")))
               (for-each (lambda (file)
                           (install-file file include))
                         (scandir "." (lambda (name) (string-match "\\.h$" name))))
               #t))))))
    (native-inputs `(("pkg-config" ,pkg-config)))
    (inputs
     `(("htslib" ,my/htslib)
       ("ncurses" ,ncurses)
       ("perl" ,perl)
       ("python" ,python)
       ("zlib" ,zlib)))
    (home-page "http://samtools.sourceforge.net")
    (synopsis "Utilities to efficiently manipulate nucleotide sequence alignments")
    (description
     "Samtools implements various utilities for post-processing nucleotide
sequence alignments in the SAM, BAM, and CRAM formats, including indexing,
variant calling (in conjunction with bcftools), and a simple alignment
viewer.")
    (license expat)))

(define %pigx-bsseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-bsseq
  (package
    (name "pigx_bsseq")
    (version %pigx-bsseq-version)
    (source (string-append (getcwd) "/pigx_bsseq-" version ".tar.gz"))
    (build-system gnu-build-system)
    (arguments
     `(#:phases
       (modify-phases %standard-phases
         (add-before 'check 'set-timezone
           ;; The readr package is picky about timezones.
           (lambda* (#:key inputs #:allow-other-keys)
             (setenv "TZ" "UTC+1")
             (setenv "TZDIR"
                     (string-append (assoc-ref inputs "tzdata")
                                    "/share/zoneinfo"))
             #t))
         (add-after 'install 'wrap-executable
           ;; Make sure the executable finds all R modules.
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (wrap-program (string-append out "/bin/pigx-bsseq")
                 `("R_LIBS_SITE" ":" = (,(getenv "R_LIBS_SITE")))
                 `("PYTHONPATH"  ":" = (,(getenv "PYTHONPATH")))))
             #t)))))
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)
       ("tzdata" ,tzdata)))
    (inputs
     `(("r-minimal" ,r-minimal)
       ("r-annotationhub" ,r-annotationhub)
       ("r-dt" ,r-dt)
       ("r-genomation" ,r-genomation)
       ("r-methylkit" ,r-methylkit)
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
       ("samtools" ,my/samtools)))
    (home-page "https://github.com/BIMSBbioinfo/pigx_bsseq/")
    (synopsis "Bisulfite sequencing pipeline from fastq to methylation reports")
    (description "PiGx BSseq is a data processing pipeline for raw fastq
read data of bisulfite experiments; it produces reports on aggregate
methylation and coverage and can be used to produce information on
differential methylation and segmentation.")
    (license gpl3+)))

pigx-bsseq
