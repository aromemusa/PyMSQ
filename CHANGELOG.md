# Changelog

## [0.1.2] - 2026-02-14
### Fixed
- Corrected internal genotype recoding in parent-specific marker-effect computation: the previous implementation could produce incorrect -1/0/1 codes when the major allele was coded as 0.
- Revised `calculate_par_spec_me` to fix major-allele handling for robust parent-specific marker-effect encoding.
- Revised `calculate_mspar_spec_me` to improve robustness of parent-specific Mendelian sampling marker-effect calculations.
- Revised `trait_matrices` to force marker-effect columns to be contiguous so that `msvarcov` no longer fails for single-trait cases.

### Documentation
- Updated README/tutorial citation guidance to cite the PyMSQ software paper (BMC Bioinformatics, 2026; DOI: 10.1186/s12859-026-06392-5) and the similarity-matrix method paper (JABG, 2025; DOI: 10.1111/jbg.12930).
- Added funding statement.