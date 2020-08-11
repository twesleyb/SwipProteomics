# Manuscript Outline

* 00_TITLE.tex
* 01_ABSTRACT.tex
* 02_INTRO.tex
* 03_RESULTS.tex - Contains the results as well as figures.
* 04_DISCUSSION.tex
* 05_ACKNOWLEDGEMENTS.tex
* 06_CONTRIBUTIONS.tex
* 07_AUTHOR_INTERESTS.tex
* 08_RESOURCE_AVAILABILITY.tex
* 09_METHODS.tex
* 10_REFERENCES.bib
* 11_FUNDING_SOURCES.tex

Compile the main tex document with `pdflatex` or `pdflatex-quiet`.
```
pdflatex manuscript.tex
```
Produces: `manuscript.pdf`.

#### Other important files:
* `elife.cls` Specifies the formatting of the document.
* `vancouver-elife.bst` Specifies the formatting of the BibTeX `.bib` file.
