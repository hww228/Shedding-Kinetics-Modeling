# Exponential Decay Model with ct_value Transformation

## Introduction

This project demonstrates a Bayesian modeling approach using Stan to estimate the exponential decay of RNA concentration based on qPCR cycle threshold (ct_value) data. In this model, ct_value measurements are first transformed into RNA concentrations via a linear transformation before fitting an exponential decay model.

## Model Overview

The model assumes that the RNA concentration on the log$_{10}$ scale decays exponentially over time. Because the data are collected as ct_value measurements, we first convert these values into RNA concentrations using a linear model with parameters beta0 and beta1. The mathematical formulation is as follows:


## Repository Structure
-   data/ : Contains data **wang2020fecal**. A data.xlsx → Corresponds to ORF1ab gene measurements; B data.xlsx → Corresponds to N gene measurements

-   code/reproducibility.R: Contains code to reproduce case 1: Single Subject without Censored Data 
-   code/test.R: Contains code for testing data **kim2020viral**: Single Subject without Censored ct_value. 
-   code/SinglewoCensor_wang2020fecal.R: Contains code for testing data **wang2020fecal**: Treated as Single Subject without Censored ct_value. 
-   code/CtViralcombo.R: Contains code for testing data **wang2020fecal** and **wolfel2020virological**: Treated as Multiple Subjects without Censored ct_value and viral load data. 
-   models/expon1ct.stan: Contains the Stan model with ct_value.
-   models/expon1_combined.stan: Contains the Stan model with both ct_value and viral load data.

## Prerequisites

Make sure you have R and the following packages installed first:

``` r
install.packages(c("dplyr", "readxl", "rstan"))
```

## Using pre-duilt Docker image

We have built an image with Rstudio (with R version 4.4.3) and rstan. Users can use the following docker pull command:

```
docker pull ywan446/shedding-hub
```

and run the container with

```
docker run -e PASSWORD="<YOURPASSWORD>" -p 8787:8787 ywan446/shedding-hub:1.0
```

Then, users can go to "http://localhost:8787/" and type in 

```
USERNAME: rstudio
PASSWORD: <YOURPASSWORD>
```