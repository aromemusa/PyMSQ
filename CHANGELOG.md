# Changelog

## [0.1.3] - 2026-02-24
### Fixed
- Corrected `calculate_par_spec_me` genotype recoding to use a consistent major/minor convention and to handle both haplotype input (2× rows per individual) and common genotype codings (e.g., 0/1/2 and −1/0/+1) safely.
- Replaced deprecated `pkg_resources` data loading with `importlib.resources` for packaged data access.