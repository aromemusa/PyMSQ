# PyMSQ: documentation

PyMSQ is a python package for estimating Mendelian sampling-related quantities. PyMSQ is composed of one module [`msq`](docs/documentation_msq.md) for simplicity. 

Include the following code in your python file to use its functions:

`from PyMSQ import msq`

PyMSQ can be installed using `pip` with the command

`pip install PyMSQ`

This documentation describes the main functions of the PyMSQ package:
* Derivation of population-specific covariance matrices reflecting within-family linkage disequilibrium
* Estimation of Mendelian sampling variance for traits and aggregate genotype for multiple traits, as well as their covariances
* Derivation of similarity matrices based on Mendelian sampling values between individuals
* Estimation of selection criteria for genetic evaluation of plants and animals

A tutorial on how to use PyMSQ's functions can be found [`here`](docs/Illustration of PyMSQ functions.md).
