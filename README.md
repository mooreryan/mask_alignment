# Mask Alignment

Mask gap columns in sequence alignments.

Decently fast, so you can use it with big alignment files.

## Install

If you are an OCaml developer, clone the repo, then run `opam install .`.

If you are not, instructions to set up an OCaml development environment can be found [here](https://ocaml.org/learn/tutorials/up_and_running.html) or [here](https://dev.realworldocaml.org/install.html).

### Run tests

If you want, you can run the tests. From the root of the repository, run:

```
$ dune runtest
```

## Usage

```
$ mask_alignment <aln.fa> <max_gap_percent> > <output.fa>
```

## Example

Take an alignment file `silly.aln.fa`, and mask (a.k.a., remove) any column with >= 95% gaps.

```
$ mask_alignment silly.aln.fa 95 > silly.aln.masked.fa
```

## License

[![license MIT or Apache
2.0](https://img.shields.io/badge/license-MIT%20or%20Apache%202.0-blue)](https://github.com/mooreryan/pasv)

Copyright (c) 2021 - 2023 Ryan Moore.

Licensed under the Apache License, Version 2.0 or the MIT license, at your option. This program may not be copied, modified, or distributed except according to those terms.
