# wieczorek-thesis-code
Simulations and data analyses from Jerzy Wieczorek's PhD thesis:

Wieczorek, J. (2018) *Model Selection and Stopping Rules for High-Dimensional Forward Selection.* Doctoral dissertation. Carnegie Mellon University, Department of Statistics & Data Science.

## Main simulations
The `MainSims` folder contains code (and saved simulation results) to replicate the main simulations from Chapter 6, including Figures 6.1 through 6.5.

These simulations rely on functions defined (and libraries loaded) in `Functions/CV_Functions.R`.

## Examples using real data
The `Examples` folder contains code and data to replicate the real-data examples from Chapter 6, including Figures 6.6 and 6.7 and Table 6.1.

We include a copy of the prostate dataset from the `Data` tab of the [**Elements of Statistical Learning** website](https://web.stanford.edu/~hastie/ElemStatLearn/).

The Boston housing dataset is not included directly here but is loaded through the `MASS` package in `R`.

For the Million Song Dataset, the 201MB data file is not included, but `MSD_EDA.R` contains instructions to download and pre-process it.

## Other smaller simulations
The `SmallSims` folder contains code (and saved simulation results) to replicate other simulations from throughout the thesis: Figures 4.1, 4.2, 5.1, 5.2, and 7.1.

## `R` packages used for the thesis
The primary packages used are `leaps` for implementing Forward Selection, `doParallel` for parallel computing of the main simulations, and `ggplot2` for generating most of the figures.

A complete list of the `R` packages called in these scripts is: `doParallel`, `ggplot2`, `grid`, `gridExtra`, `Hmisc`, `leaps`, `MASS`, `Matrix`, `matrixStats`, `plotrix`, `plyr`, `R.utils`, and `RColorBrewer`.
