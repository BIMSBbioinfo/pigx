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
;;; v1.0.1-291-g164cb4da63

(use-modules (guix packages) (gnu packages))

(define %pigx-rnaseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-rnaseq-development
  (package (inherit (specification->package "pigx-rnaseq"))
    (version %pigx-rnaseq-version)
    (source (string-append (getcwd) "/pigx_rnaseq-" version ".tar.gz"))
    (inputs
     `(("hisat2" ,(specification->package "hisat2"))
       ("fastp" ,(specification->package "fastp"))
       ("python-deeptools" ,(specification->package "python-deeptools"))
       ("gprofiler2", (specification->package "r-gprofiler2"))
       ,@(package-inputs (specification->package "pigx-rnaseq"))))
    (native-inputs
     `(("autoconf" ,(specification->package "autoconf"))
       ("automake" ,(specification->package "automake"))))))

pigx-rnaseq-development
