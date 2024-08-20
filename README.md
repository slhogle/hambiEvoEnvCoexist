# hambiBottomUp

[Click here to view rendered notebooks of the analysis.](https://hambiloopairwisenmr-slhogl-c8d6c219bb8d719ac159d67b545bcea32595.utugit.fi/)

## Manuscript:

"TBD"

[Preprint available from bioRxiv]()

## Data

Data and code here is provided under GPL3. Feel free to use or remix as you see fit.

[![DOI](https://zenodo.org/badge/.svg)](https://zenodo.org/badge/latestdoi/)

### Project structure
1.  `/R` contains R scripts
2.  `/data` contains data that has been processed in some way for later use
3.  `/_data_raw` contains unprocessed data scraped from compute cluster
4.  `/figs` contains figures generated from R scripts

## Availability

The rendered project site is available at <https://hambiloopairwisenmr-slhogl-c8d6c219bb8d719ac159d67b545bcea32595.utugit.fi/>. The website has been produced using [Quarto notebooks](https://quarto.org/).

This GitLab repository (<https://gitlab.utu.fi/slhogl/hambiLOOPairwiseNMR>) hosts the code and data for this project. The rendered webpage can be fully recreated using the code at <https://gitlab.utu.fi/slhogl/hambiLOOPairwiseNMR>. 

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. To get up and running you can do:

``` r
install.packages("renv")
renv::restore()
```
