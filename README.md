# Bayesian Spatial Modeling with Gaussian Processes

## Overview

This project focuses on the spatial analysis of air pollution across European locations using a Bayesian Gaussian Process (GP) model.

The objective is to investigate spatial dependencies in PM10 concentrations, estimate the underlying spatial process, and perform spatial predictions while accounting for uncertainty.

The project was developed as part of a university assignment at Politecnico di Torino.

## Methodology

A hierarchical Bayesian spatial model was implemented to describe the logarithm of PM10 concentrations as a function of geographical coordinates and a spatially correlated latent process.

The model incorporates:

- A linear spatial trend based on geographical coordinates.
- A Gaussian Process with an exponential covariance function.
- A Gaussian observation error term.
- Weakly informative prior distributions for model parameters.

Posterior inference was performed using a custom Markov Chain Monte Carlo (MCMC) algorithm combining:

- **Gibbs Sampling** for parameters with analytically available full conditional distributions.
- **Adaptive Metropolis-Hastings** for the spatial correlation parameter.

## Analysis and Results

The project includes:

1. **Exploratory Spatial Data Analysis:** visualization of PM10 concentrations across European locations to investigate spatial patterns.

2. **Bayesian Parameter Estimation:** estimation of regression coefficients, variance components, and spatial correlation parameters.

3. **MCMC Diagnostics:** evaluation of parameter traces and autocorrelation to assess sampling behavior.

4. **Spatial Dependence Analysis:** estimation of the covariance function and its 95% posterior credible intervals.

5. **Spatial Prediction:** estimation of the latent Gaussian Process over a regular geographical grid, including posterior means and credible intervals.

## Technologies

- R
- Bayesian Statistics
- Gaussian Processes
- Markov Chain Monte Carlo (MCMC)
- Spatial Data Analysis
- ggplot2

MSc in Mathematical Engineering

Politecnico di Torino
