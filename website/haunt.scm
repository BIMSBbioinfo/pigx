;; -*- geiser-scheme-implementation: guile -*-

(use-modules (haunt builder blog)
             (haunt builder assets)
             (haunt reader)
             (haunt reader commonmark)
             (haunt site)
             (haunt post)
             (haunt page)
             (haunt html)
             (srfi srfi-1))


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
               ,body
               (footer
                (img (@ (alt "decoration")
                        (class "end-decoration")
                        (src "images/end.svg")))
                (p "Made by the Akalin lab at the Berlin Institute of Medical Systems Biology in 2018"))))))

(define default-theme
  (theme #:name "Default"
         #:layout default-layout
         #:post-template
         (lambda (post)
           `((h1 ,(post-ref post 'title))
             ,(post-sxml post)))))


;;; XXX This is needed because Guile Commonmark does not yet support
;;; raw HTML blocks.  See
;;; https://github.com/OrangeShark/guile-commonmark/issues/8
(use-modules (commonmark)
             (haunt post)
             (haunt reader)
             (srfi srfi-1)
             (sxml transform))

;; Render markdown but replace any "###### INSERT EXAMPLES HERE" with
;; the examples.
(define commonmark-reader-with-examples
  (make-reader
   (make-file-extension-matcher "md")
   (lambda (file)
     (call-with-input-file file
       (lambda (port)
         (values (read-metadata-headers port)
                 (pre-post-order
                  (commonmark->sxml port)
                  `((h6 . ,(lambda (sym tree)
                             (if (and (string? tree)
                                      (string=? "INSERT EXAMPLES HERE" tree))
                                 '(div (@ (id "example-scroller"))
                                       (div (@ (id "examples"))
                                            (img (@ (class "example")
                                                    (src "images/report1.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report2.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq1.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq2.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq3.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq4.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq5.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq6.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq7.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq8.png")
                                                    (alt "Example report screenshot")))
                                            (img (@ (class "example")
                                                    (src "images/report-rnaseq9.png")
                                                    (alt "Example report screenshot")))))
                                 tree)))
                    (*text* . ,(lambda (sym text) text))
                    (*default* . ,(lambda arg arg))))))))))

(define* (alias-post file-name alias #:key theme)
  "Return a builder procedure that builds FILE-NAME as ALIAS."
  (lambda (site posts)
    (make-page alias
               (render-post theme site
                            (find (lambda (post)
                                    (string=? (post-file-name post) file-name))
                                  posts))
               sxml->html)))

(site #:title "PiGx: Pipelines in Genomics"
      #:domain "http://bioinformatics.mdc-berlin.de/pigx"
      #:default-metadata
      '((author . "BIMSB Bioinformatics")
        (email  . "ricardo.wurmus@mdc-berlin.de"))
      #:readers (list commonmark-reader-with-examples html-reader)
      #:builders (list (blog #:theme default-theme #:collections '())
                       (alias-post "posts/index.md" "index.html" #:theme default-theme)
                       (static-directory "static" ".")))
