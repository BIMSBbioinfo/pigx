;;; PiGx_scrnaseq - reports pipeline for single cell RNAseq experiments.
;;; Copyright © 2018, 2019 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx_scrnaseq.
;;;
;;; PiGx_scrnaseq is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PiGx_scrnaseq:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; To install the package from a release tarball do this:
;;;
;;;  $ guix package --with-source=pigx_scrnaseq-0.0.1.tar.gz -f guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.14.0-3089-g84c195e50

(use-modules (guix packages) (gnu packages))

(define %pigx-scrnaseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-scrnaseq-development
  (package (inherit (specification->package "pigx-scrnaseq"))
    (version %pigx-scrnaseq-version)
    (source (string-append (getcwd) "/pigx_scrnaseq-" version ".tar.gz"))
    (native-inputs
     `(("autoconf" ,(specification->package "autoconf"))
       ("automake" ,(specification->package "automake"))))))

pigx-scrnaseq-development
