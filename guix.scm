;;; PiGx_scrnaseq - reports pipeline for single cell RNAseq experiments.
;;; Copyright © 2018 Ricardo Wurmus <rekado@elephly.net>
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
;;; v0.13.0-4757-ga0837294a

(use-modules (guix packages)
             (guix licenses)
             (guix download)
             (guix build-system ant)
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

;; FIXME: This package includes pre-built jars.
(define-public dropseq-tools
  (package
    (name "dropseq-tools")
    (version "1.13")
    (source
     (origin
       (method url-fetch)
       (uri "http://mccarrolllab.com/download/1276/")
       (file-name (string-append "dropseq-tools-" version ".zip"))
       (sha256
        (base32
         "0yrffckxqk5l8b5xb6z4laq157zd9mdypr2p4b4vq2bhjzi1sj0s"))
       ;; Delete bundled libraries
       (modules '((guix build utils)))
       (snippet
       	'(begin
	   ;; TODO: not all libraries are available in Guix yet
	   ;;(for-each delete-file (find-files "jar/lib" "\\.jar$"))
	   (delete-file-recursively "3rdParty")))))
    (build-system ant-build-system)
    (arguments
     `(#:tests? #f ; test data are not included
       #:test-target "test"
       #:build-target "all"
       #:source-dir "public/src/"
       #:jdk ,icedtea-8
       #:make-flags
       (list (string-append "-Dpicard.executable.dir="
			    (assoc-ref %build-inputs "java-picard") "/share/java/"))
       #:phases
       (modify-phases %standard-phases
         ;; There is no installation target
         (replace 'install
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let* ((out     (assoc-ref outputs "out"))
                    (bin     (string-append out "/bin"))
                    (share   (string-append out "/share/dropseq-tools/"))
                    (scripts (list "CreateIntervalsFiles"
				   "ConvertToRefFlat"
				   "GatherGeneGCLength"
				   "GatherReadQualityMetrics")))
               (for-each mkdir-p (list bin share))
	       (install-file "dist/dropseq.jar" share)
	       (for-each (lambda (script)
			   (chmod script #o555)
			   (install-file script bin))
			 scripts)
	       (substitute* (map (lambda (script)
				   (string-append bin "/" script))
				 scripts)
		 (("^java") (which "java"))
		 (("jar_deploy_dir=.*")
		  (string-append "jar_deploy_dir=" share "\n")))))))))
    (inputs
     `(("java" ,icedtea-8)
       ("java-picard" ,java-picard)
       ("java-log4j-1.2" ,java-log4j-1.2-api)
       ("java-commons-math3" ,java-commons-math3)
       ("java-commons-jexl2" ,java-commons-jexl-2)
       ("java-commons-collections4" ,java-commons-collections4)
       ("java-commons-lang2" ,java-commons-lang)
       ("java-commons-io" ,java-commons-io)
       ("java-snappy-1.0.3-rc3" ,java-snappy-1)
       ("java-slf4j-api" ,java-slf4j-api)
       ("java-guava" ,java-guava) ; TODO: needs version 22, but we only have 20
       ;; la4j 0.5.5
       ;; biojava-core 4.0.0
       ;; biojava-alignment 4.0.0
       ;; jdistlib
       ;; slf4j-log4j12-1.7.22
       ;; simple-xml-2.0.4
       ("java-snakeyaml" ,java-snakeyaml)))
    (native-inputs
     `(("unzip" ,unzip)
       ("java-testng" ,java-testng)))
    (home-page "http://mccarrolllab.com/dropseq/")
    (synopsis "Tools for Drop-seq analyses")
    (description "Drop-seq is a technology to enable biologists to
analyze RNA expression genome-wide in thousands of individual cells at
once.  This package provides tools to perform Drop-seq analyses.")
    (license expat)))

(define %pigx-scrnaseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-scrnaseq
  (package
    (name "pigx_scrnaseq")
    (version %pigx-scrnaseq-version)
    (source (string-append (getcwd) "/pigx_scrnaseq-" version ".tar.gz"))
    (build-system gnu-build-system)
    (arguments
     `(#:tests? #f ; requires network access
       #:phases
       (modify-phases %standard-phases
         (add-after 'install 'wrap-executable
           ;; Make sure the executable finds all R modules.
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (wrap-program (string-append out "/bin/pigx-scrnaseq")
                 `("R_LIBS_SITE" ":" = (,(getenv "R_LIBS_SITE")))
                 `("PYTHONPATH"  ":" = (,(getenv "PYTHONPATH")))))
             #t)))))
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)))
    (inputs
     `(("dropseq-tools" ,dropseq-tools)
       ("python-wrapper" ,python-wrapper)
       ("python-pyyaml" ,python-pyyaml)
       ("python-pandas" ,python-pandas)
       ("python-numpy" ,python-numpy)
       ("ghc-pandoc" ,ghc-pandoc)
       ("ghc-pandoc-citeproc" ,ghc-pandoc-citeproc)
       ("snakemake" ,snakemake)
))
    (home-page "https://github.com/BIMSBbioinfo/pigx_scrnaseq/")
    (synopsis "TODO")
    (description "TODO")
    (license gpl3+)))

pigx-scrnaseq
