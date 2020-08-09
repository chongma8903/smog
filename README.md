# smog
[![Travis-CI Build Status](https://travis-ci.org/chongma1989/smog.svg?branch=master)](https://travis-ci.org/chongma1989/smog)

Structural Modeling by using Overlapped Group Penalty

## Installation
```r
# Install smog from CRAN
install.packages("smog")

# or install the source type package from GitHub:
# install.packages("devtools")
devtools::install_github("chongma8903/smog")
```

## Features
* fits a linear non-penalized phenotype (demographic) variables and penalized groups of prognostic effect and predictive effect.
* satisfies such hierarchy structures that if a predictive effect exists, its prognostic effect must also exist.
* can deal with continuous, binomial or multinomial, and survival response variables.
* incorporates the iterative shrinkage-thresholding algorithm (ISTA) and the alternating direction method of multipliers algorithms (ADMM).

## Usage
Create a new S3 class of `smog`, and the kernal function `smog.default` (or `smog.formula`) returns an object of the S3 class `smog`. 

* `smog.default` input the data and parameters to yield a model of the class `smog`.
* `smog.formula` can accept `formula` to fit the model for the data.
* `predict.smog` produces the predicted response values for new data, provided a fitted model. 
* `cv.smog` provides cross-validation analysis based on the data.  
* `plot.smog` displays a panel of three plots to reflect the convergence performance in the algorithm.
* `print.smog` outputs the coefficients table for the fitted model.
