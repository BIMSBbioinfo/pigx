;;; PiGx SARS-CoV2 wastewater sequencing pipeline
;;; Copyright © 2021 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx SARS-CoV2 wastewater sequencing pipeline
;;;
;;; This is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment:
;;;
;;;  $ guix environment -m manifest.scm
;;;
;;; When the shell variable USE_GUIX_INFERIOR is set (to any value), a
;;; specific set of Guix channels will be used to build the
;;; environment.

(use-modules
 (guix profiles)
 (guix channels)
 (guix inferior)
 (ice-9 match))

(define channels
  (list (channel
         (name 'guix-bimsb)
         (url "https://github.com/BIMSBbioinfo/guix-bimsb.git")
         (commit
          "be42b099867a29a68a9ffbf304dc7f9137b75fe3"))
        (channel
         (name 'guix)
         (url "https://git.savannah.gnu.org/git/guix.git")
         (commit
          "01e33a031e493477d930b9383d397fea012a3b1a")
         (introduction
          (make-channel-introduction
           "9edb3f66fd807b096b48283debdcddccfea34bad"
           (openpgp-fingerprint
            "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA"))))))

(define (lookup name)
  (specification->package name))

(define (lookup-inferior name)
  (define inferior
    (inferior-for-channels channels))
  (match (lookup-inferior-packages inferior name)
    ((first . rest) first)
    (_ (error
        (format #false "Could not find package `~a'.~%" name)))))

(define %packages
  (list "bwa"
        "ensembl-vep"
        "fastqc"
        "kraken2"
        "krona-tools"
        "lofreq"
        "prinseq"
        "samtools"
        "snakemake"
        "r-minimal"
        "r-plotly"
        "r-rmarkdown"
        "python-wrapper"
        "python-pyyaml"))

(define %native-packages
  (list "autoconf"
        "automake"))

(packages->manifest
 (let ((how (if (getenv "USE_GUIX_INFERIOR")
                lookup-inferior lookup)))
   (map how (append %packages %native-packages))))
