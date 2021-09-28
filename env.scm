;; For a persistent profile with development packages use
;;    guix package -p dev -m env.scm
(use-modules (guix packages))
(packages->manifest
 (append (map cadr (package-direct-inputs
                    (specification->package "pigx-chipseq")))
         (map specification->package '("autoconf" "automake"))
		 ;; additional dependencies
         (map specification->package '("r-rsubread" "samblaster"))
         ;; add pacages for Deseq Report
         (map specification->package '(
                                       "r-ggrepel"
                                       "r-deseq2"
                                       "r-dt"
                                       "r-pheatmap"
                                       "r-corrplot"
                                       "r-reshape2"
                                       "r-scales"
                                       "r-crosstalk"
                                       "r-gprofiler2"
                                       "r-summarizedexperiment"
                                       "r-hexbin"
                                       ))
		 ))
