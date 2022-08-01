# CHVfilter (Compressed hierarchical Vecchia filter)

Reference: Chakraborty, A., & Katzfuss, M. (2022). High-dimensional nonlinear spatio-temporal filtering by compressing hierarchical sparse Cholesky factors.

We propose a scalable filtering approach based on a hierarchical sparse Cholesky representation of the filtering covariance matrix. At each time point, we  compress the sparse Cholesky factor into a dense matrix with a small number of columns. After applying the evolution to each of these columns, we decompress to obtain a hierarchical sparse Cholesky factor of the forecast covariance, which can then be updated based on newly available data.

Note - This repository contains a directory named GPvecchia which exactly resembles the R package `GPvecchia`. However, this local directory contains few additional codes needed for running our filtering algorithm. Hence, we must load this directory as well (this will be automatically implemented if you follow the instructions below).

## Installation

At first, clone this directory from git in your local library using the following command

```
git clone https://github.com/katzfuss-group/CHVfilter.git
```

## Usage

### 1\. Installing prerequisite packages

After the repository is cloned, open the `install-all-packages.R` file in your local repository and run all the lines to install the prerequisite packages.

### 2\. Generating plots that do not require simulations

The codes for the plots that do not need simulations are stored in the file `no-simulation-plots.R`. The R code is sectioned properly to indicate several plots. You may run the whole file to generate all the plots, or run sections separately to generate selective plots.

### 3\. Running the comparison with EKVL filter

To run the log-score comparison with EKVL filter, go to the `filtering-lorenz-chv-ekvl` directory. This directory contains two `R` files. At first, open the `run-filtering-lorenz-chv-ekvl.R` file, and run all the lines to generate the simulations from both EKVL and CHV filtering. This should take several hours to complete based on the system, often approximately four to five hours. 

Then, open the `plots-lorenz-ekvl-ehv.R` file, and run all the lines. This will generate the plots for comparing log score of CHV filter with that of EKVL filter.

### 3\. Running the comparison with RRUKF

To run the log-score comparison with EKVL filter, go to the `filtering-lorenz-chv-rrukf` directory. This directory contains two `R` files. At first, open the `run-filtering-lorenz-chv-redrank.R` file, and run all the lines to generate the simulations from both EKVL and CHV filtering. This should take several hours to complete based on the system, often approximately four to five hours.

Then, open the `plots-lorenz-redrank-chv.R` file, and run all the lines. This will generate the plots for comparing log score of CHV filter with that of EKVL filter.
