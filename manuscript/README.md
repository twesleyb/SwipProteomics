# Manuscript Outline

* `0_TITLE.tex` The manuscript's TITLE.
* `1_ABSTRACT.tex` The manuscript's ABSTRACT.
* `2_INTRO.tex` The manuscript's INTRODUCTION.
* `3_RESULTS.tex` The manuscript's RESULTS.
* `4_DISCUSSION.tex` The manuscript's DISCUSSION.
* `5_METHODS.tex` The manuscript's METHODS.
* `6_ACKNOWLEDGMENTS.tex` The manuscript's ACKNOWLEDGEMENTS.
* `7_FUNDING.tex` Funding sources for the work.
* `8_REFS.bib` References in bibtex format.
* `ARCHIVE/` Directory for old files.

Compile the main tex document with `pdflatex` or `pdflatex-quiet`.
```
pdflatex MAIN.tex
```
Produces: `MAIN.pdf`.

#### Other important files.
* `elife.cls` Specifies the formatting of the document.
* `vancouver-elife.bst` Specifies the formatting of the BibTeX `.bib` file.
