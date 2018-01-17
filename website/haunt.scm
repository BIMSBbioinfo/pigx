;; -*- geiser-scheme-implementation: guile -*-

(use-modules (haunt builder blog)
             (haunt builder assets)
             (haunt reader)
             (haunt reader commonmark)
             (haunt site)
             (haunt post))


(define (default-layout site title body)
  `((doctype "html")
    (head
     (meta (@ (http-equiv "Content-Type") (content "text/html; charset=UTF-8")))
     (meta (@ (http-equiv "Content-Language") (content "en")))
     (meta (@ (name "author") (content "BIMSB Bioinformatics")))
     (meta (@ (name "viewport") (content "width=device-width")))
     (title ,(or title (if site (site-title site) "PiGx: Pipelines in Genomics")))
     (link (@ (rel "stylesheet")
              (media "screen")
              (type "text/css")
              (href "css/screen.css")))
     (link (@ (rel "shortcut icon")
              (href "https://bioinformatics.mdc-berlin.de/pigx/favicon.ico"))))
    (body (@ (id "top"))
          (div (@ (id "index") (align "center"))
               (img (@ (alt "logo")
                       (src "images/logo.svg")
                       (width "30%")
                       (height "30%"))))
          (div (@ (id "page"))
               ,body))))

(define default-theme
  (theme #:name "Default"
         #:layout default-layout
         #:post-template
         (lambda (post)
           `((h1 ,(post-ref post 'title))
             ,(post-sxml post)))))

(site #:title "PiGx: Pipelines in Genomics"
      #:domain "http://bioinformatics.mdc-berlin.de/pigx"
      #:default-metadata
      '((author . "BIMSB Bioinformatics")
        (email  . "ricardo.wurmus@mdc-berlin.de"))
      #:readers (list commonmark-reader html-reader)
      #:builders (list (blog #:theme default-theme #:collections '())
                       (static-directory "static" ".")))
