# Exponential Decay Model with ct_value Transformation

## Introduction

This project demonstrates a Bayesian modeling approach using Stan to estimate the exponential decay of RNA concentration based on qPCR cycle threshold (ct_value) data. In this model, ct_value measurements are first transformed into RNA concentrations via a linear transformation before fitting an exponential decay model.

## Model Overview

The model assumes that the RNA concentration on the log$_{10}$ scale decays exponentially over time. Because the data are collected as ct_value measurements, we first convert these values into RNA concentrations using a linear model with parameters beta0 and beta1. The mathematical formulation is as follows:


## Repository Structure
-   code/reproducibility.R: Contains code to reproduce case 1: Single Subject without Censored Data 
-   code/test.R: Contains code for testing data **kim2020viral**: Single Subject without Censored ct_value. 
-   code/test2.R: Contains code for testing data **wang2020fecal**: Treated as Single Subject without Censored ct_value. 
-   models/expon1ct.stan: Contains the Stan model with ct_value.

## Prerequisites

Make sure you have R and the following packages installed first:

``` r
install.packages(c("dplyr", "readxl", "rstan"))
```
