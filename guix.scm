;;; PiGx SARS-CoV2 wastewater sequencing pipeline
;;; Copyright © 2021 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx SARS-CoV2 wastewater sequencing pipeline
;;;
;;; This is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; This environment file was developed for Guix with this list of
;;; channels:
;;;
;;; (list (channel
;;;        (name 'guix)
;;;        (url "https://git.savannah.gnu.org/git/guix.git")
;;;        (commit
;;;         "e2ab6fb0dd3e2bcaadc0d349c96c3689db5ac670")
;;;        (introduction
;;;         (make-channel-introduction
;;;          "9edb3f66fd807b096b48283debdcddccfea34bad"
;;;          (openpgp-fingerprint
;;;           "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))

(use-modules
 (guix build-system gnu)
 (guix packages)
 (guix licenses)
 (gnu packages))

(define %version
  (symbol->string (with-input-from-file "VERSION" read)))

(define (p name)
  (let ((pkg (specification->package name)))
    `(,name ,pkg)))

(define-public pigx-sars-cov2-ww-development
  (package
    (name "pigx-sars-cov2-ww")
    (version %version)
    (source
     (string-append (getcwd) "/pigx_sars-cov2-ww-" version ".tar.gz"))
    (build-system gnu-build-system)
    (inputs
     (map p (list "fastqc"
                  "kraken2"
                  "prinseq"
                  "snakemake"
                  "r-minimal"
                  "r-plotly"
                  "python-wrapper"
                  "python-pyyaml")))
    (native-inputs
     (map p (list "autoconf"
                  "automake")))
    (home-page "http://bioinformatics.mdc-berlin.de/pigx/")
    (synopsis "SARS-CoV2 wastewater sequencing pipeline")
    (description "")
    (license gpl3+)))

pigx-sars-cov2-ww-development
