## Supplemental Methods

This is very complicated.
The main file is `supplement.Rnw`. This file is compiled with `knitr` to
generate `supplement.tex`. This file must then be parsed by `latex` and
`bibtex`. In order for references and citations to work properly, you must
compile the `.tex` document twice. The order of operations is:

* knitr 
* pdflatex
* bibtex
* pdflatex
* pdflatex

The main `.Rnw` file calls R source code in `supplement.R`.
The main `.Rnw` file sources figures from `figures.tex`.
The main `.Rnw` file sources references from `bibliography.bib`.

For simplicity, all the figures and their captions are managed in a seperate
file, `figures.tex`. The raw `*.pdf` files are in `figs/`.
