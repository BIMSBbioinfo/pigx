;; For a persistent profile with development packages use
;;    guix package -p dev -m env.scm
(use-modules (guix packages))
(packages->manifest
 (append (map cadr (package-direct-inputs
                    (specification->package "pigx-chipseq")))
         (map specification->package '("autoconf" "automake"))))
