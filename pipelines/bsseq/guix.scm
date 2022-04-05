;;; PiGx_bsseq - reports pipeline for reads from bisulfite experiments.
;;; Copyright © 2017, 2018, 2020 Ricardo Wurmus <rekado@elephly.net>
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
;;; 1614862

(use-modules (guix packages) (gnu packages))

(define %pigx-bsseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-bsseq-development
  (package (inherit (specification->package "pigx-bsseq"))
    (version %pigx-bsseq-version)
    (source (string-append (getcwd) "/pigx_bsseq-" version ".tar.gz"))
    (native-inputs
     `(("autoconf" ,(specification->package "autoconf"))
       ("automake" ,(specification->package "automake"))))))

pigx-bsseq-development
