# PyMSQ: Estimating Mendelian Sampling-Related Quantities in Python
PyMSQ is an open-source Python package that enables breeders, geneticists, and quantitative biologists to estimate Mendelian sampling–related metrics—including variance, covariance, and haplotype-based similarities—in both plant and animal species. For simplicity, PyMSQ consists of a single module, [`msq`](docs/documentation_msq.md). 

## Key Features
- Within-Family Covariance
Constructs population-specific covariance matrices that capture within-family linkage disequilibrium, reflecting recombination patterns and phased marker data.

- Mendelian Sampling (Co)Variance
Estimates Mendelian sampling variance (MSV) for single or multiple traits, as well as covariances (MSCs), crucial for maintaining genetic diversity and controlling inbreeding.

- Similarity Matrices
Computes haplotype-based similarity matrices between individuals (or zygotes), focusing on shared heterozygous segments that drive within-family genetic variation.

- Selection Criteria
Offers functions to derive selection strategies (e.g., GEBVs, usefulness criteria, index-based approaches) that leverage MSV/MSC or similarity measures.

## Installation
PyMSQ is available on PyPI and can be installed via:


```python
pip install PyMSQ
```
    

## Basic Usage
Below is a minimal example illustrating how to import PyMSQ and call its core functions:


```python
from PyMSQ import msq  # Imports the msq module

# Example: Loading an included dataset
data = msq.load_package_data()

# Deriving expected LD matrices
ld_matrices = msq.expldmat(data['chromosome_data'], data['group_data'])

# Estimating Mendelian sampling (co)variances
msv = msq.msvarcov(
    gmat      = data['genotype_data'],
    gmap      = data['chromosome_data'],
    meff      = data['marker_effect_data'],
    exp_ldmat = ld_matrices,
    group     = data['group_data']
)

# Constructing similarity matrices
similarity = msq.simmat(
    gmat      = data['genotype_data'],
    gmap      = data['chromosome_data'],
    meff      = data['marker_effect_data'],
    group     = data['group_data'],
    exp_ldmat = ld_matrices
)
```


## Tutorial
A tutorial detailing each function’s parameters, usage examples, and best practices can be found [`here`](docs/Illustration_of_PyMSQ_functions.md). This tutorial walks you through:

1. **Loading** your own data or the bundled Holstein-Friesian dataset,

2. **Building** LD matrices for each chromosome,

3. **Estimating** Mendelian sampling (co)variances,

4. **Deriving** haplotype-based similarity,

5. **Applying** selection strategies using advanced metrics.


## Getting Help
- **Issues & Feature Requests**
If you encounter bugs, have feature requests, or need additional clarification, please open an issue on the [`PyMSQ GitHub repository`](https://github.com/aromemusa/PyMSQ).

- **License**
PyMSQ is released under the MIT License, allowing both academic and commercial use.



**Happy analyzing!**
We hope PyMSQ supports your work in breeding, helping you balance short-term genetic gains with the long-term preservation of essential haplotype diversity.

