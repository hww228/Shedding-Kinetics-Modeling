# Exponential Decay Model with ct_value Transformation

## Introduction

This project demonstrates a Bayesian modeling approach using Stan to estimate the exponential decay of RNA concentration based on qPCR cycle threshold (ct_value) data. In this model, ct_value measurements are first transformed into RNA concentrations via a linear transformation before fitting an exponential decay model.

## Model Overview

The model assumes that the RNA concentration on the log$_{10}$ scale decays exponentially over time. Because the data are collected as ct_value measurements, we first convert these values into RNA concentrations using a linear model with parameters beta0 and beta1. The mathematical formulation is as follows:

1.  Exponential Decay (RNA Scale)

For each observation at time $t_i$, the mean RNA concentration is modeled as: $$\mu_{\text{RNA}}(t_i) = \frac{c_0}{\ln(10)} - \frac{a_0\,t_i}{\ln(10)}$$ where:

$c_0$ is the RNA concentration (on a log scale) at time $t = 0$

$a_0$ is the decay rate.

We assume that $$\text{RNA}_i \sim \mathcal{N}\!\Biggl(\mu_{\text{RNA}}(t_i),\,\sigma^2\Biggr)$$

2.  Transformation to ct_value:

The observed ct_value is linearly related to the RNA concentration by the transformation: $$
\text{RNA}_{\text{obs}} = \beta_0 + \beta_1\, \text{ct_value}.
$$ Solving for the ct_value, we have: $$\text{ct_value} = \frac{{\text{RNA}_{\text{obs}}}- \beta_0}{\beta_1} $$

so that: $$
\mu_{\text{ct}}(t_i) = \frac{\mu_{\text{RNA}}(t_i) - \beta_0}{\beta_1}
= \frac{\frac{c_0}{\ln(10)} - \frac{a_0\,t_i}{\ln(10)} - \beta_0}{\beta_1}.
$$

$$\sigma_{\text{ct}}(t_i) = \frac{\sigma^2}{\beta_1^2}$$

3.  Likelihood

Assuming that the transformed RNA concentration follows a normal distribution on the ct_value scale, we write: $$\text{ct_value}_i \sim \mathcal{N}\!\Biggl(\mu_{\text{ct}}(t_i),\, \sigma_{\text{ct}}(t_i)\Biggr). $$ In Stan, the standard deviation is the square root of the variance. Hence, we use: $$\sigma_{\text{ct}} = \frac{\sigma}{|\beta_1|}$$

so that:

$$\text{ct_value}_i \sim \mathcal{N}\!\Bigl(\mu_{\text{ct}}(t_i),\, \frac{\sigma}{|\beta_1|}\Bigr).$$

To keep things consistent with our expectation that $\beta_1$ is negative in theory, we use the absolute value `fabs(beta1)`(or constrain $\beta_1$ to be negative and then use $-\beta_1$) to ensure that the standard deviation is positive.

## Repository Structure

-   models/expon1ct.stan: Contains the Stan model with ct_value.

## Prerequisites

Make sure you have R and the following packages installed first:

``` r
install.packages(c("dplyr", "readxl", "rstan"))
```
