# hambiABBAssembly

[Click here to view rendered notebooks of the analysis.](https://slhogle.github.io/hambiABBAssembly/)

## Manuscript:

"TBD"

[Preprint available from bioRxiv]()

## Data

Data and code here is provided under GPL3. Feel free to use or remix as you see fit.

[![DOI](https://zenodo.org/badge/.svg)](https://zenodo.org/badge/latestdoi/)

### Project structure
1.  `/R` contains R scripts and [Quarto notebooks](https://quarto.org/)
2.  `/renv` code for [`renv`](https://rstudio.github.io/renv/index.html)
3.  `/data` contains data that has been processed in some way for later use
4.  `/_data_raw` contains unprocessed data scraped from compute cluster
5.  `/figs` contains figures generated from R scripts

## Availability

The rendered project site is available at <https://slhogle.github.io/hambiABBAssembly/>. The website has been produced using [Quarto notebooks](https://quarto.org/).

This GitHub repository (<https://github.com/slhogle/hambiABBAssembly>) hosts the code and data for this project. The rendered webpage can be fully recreated using the code.

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. To get up and running you can do:

``` r
install.packages("renv")
renv::restore()
```
