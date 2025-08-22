(lang dune 3.20)

(cram enable)

(name mask_alignment)

(generate_opam_files true)

(authors "Ryan M. Moore")

(maintainers "Ryan M. Moore")

(homepage "https://github.com/mooreryan/mask_alignment")

(documentation "https://github.com/mooreryan/mask_alignment")

(bug_reports "https://github.com/mooreryan/mask_alignment/issues")

(version 1.0.2)

(license "MIT")

(package
 (name mask_alignment)
 (synopsis "Mask gap columns in sequence alignments")
 (depends
  (ocaml
   (and
    (>= "4.14.0")))
  (bio_io
   (and
    (>= "0.6.0")))
  (base
   (and
    (>= "v0.16.0")))
  (stdio
   (and
    (>= "v0.16.0")))))
