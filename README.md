Scripts used for investigating genomic signatures of climate adaptation in Eucalyptus microcarpa

Author: Rebecca Jordan 2017

-- SNP genotyping QC --

RepCheck_similarity_bylocus.py

RepCheck4.py

-- BayeScan --

Create input file:  PEDtoBayeScan.sh / .spid

Process results:  bayescan_outliers.sh, bayescan_outliers_plus_significance.sh

-- FDIST2 (Lositan) --

Create input file:  PEDtoGenepop.sh / PEDtoGENPOP.spid

Process results:  lositan_outliers_plus_significance.sh, Lositan_qval.R

-- Hierarchical FDIST2 (Arlequin) --

Create input file:  PEDtoArlequin.sh / .spid

Process results:  arl_outliers.sh, arl_outliers_plus_significance.sh, arl_qval.R

-- Bayenv Xt X --

Create input file:  PEDtoBayenv.sh / .spid

Mean covariance matrix:  meanCovMat.py, meanCovMat.sh

Process results: XtX_summary.py, bayenv_summary_to_loci.sh

-- Bayenv2 Environmental association --

Create input file: PEDtoBayenv.sh / .spid

Mean covariance matrix:  meanCovMat.py, meanCovMat.sh

Process results: Bayenv2_results.R, bayenv_results.sh (wrapper script for bayenv_enviro_summary.py, bayenv_summary_to_loci.sh and bayenv_table_results.py)

-- Project allele frequency changes under climate change --

Future_allele_frq_change.R