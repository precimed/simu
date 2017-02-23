# Synthetic genotypes and phenotypes

SynGP is about generating synthetic data.

To model  effects of linkage disequilibrium in our synthetic data we use realistic genotypes and calculate GWAS estimates of the SNP effect sizes.
These realiztic genotypes were obtained by Simulation_pipeline.py script,
which applies HapGen2 utility to 1000 Genomes data and produce genotypes of 100K individuals, 80M SNPs per each invididual (no missing calls).
This data is quite heavy (2 TB uncompressed, 139 GB in gzip).
This data is split into chromosomes and chunks; all chunk files share the same set of individuals.

We extract several subsets of this data to speedup preliminary computations. For example,

- Folder `<NORSTORE>/SYNGP/EUR_10K_2M_merged` contain first 10K individuals and 2518104 SNPs.
  Here 2518104 is a subset of 2558411 SNPs that we commonly use in matlab code.
- Folder `<NORSTORE>/SYNGP/EUR_100K_1M_merged` contains all 100K individuals, but just 1188973 SNPs.
  These SNPs is a subset of 1190321 SNPs (matlab reference file for 1M template (w_hm3)

In the code these "subsets" are refered to as "frames".
Each frame has a specific number of SNPs and genotyped subjects.
Within each frame we can store vectors of length "snps" or "subj", for example:

* truebeta - true effect sizes generated for synthetic mixture
* mixture  - component of the mixture; 0 = null component, the rest = causal components
* truepheno - true phenotype (without environmental component)
* phenotype - observed phenotype (with environmental component
* gwaspval  - gwas p-values
* gwasbeta  - gwas estimate of the effect size
* gwassize  - number of subjects per SNPS involved in GWAS
* gwaslogp   - minus decimal logarithm of p-value
* gwaszscore - signed summary statistic from GWAS
* genofreq   - genotype frequencies (AA, Aa, aa)

In addition to these vectors or matrices each frame has meta-information:
* snps - number of SNPs in the frame
* subj - number of subjects in the frame
* name - the name of the frame (for example 'EUR_10K_2M_merged')
* opts - misc options used to produce this frame

It is possible to move data from one frame to another with 'reframe' function.
It is possible to create frames and fill them with data using these methods:

- make_empty_frame
- make_truebeta_gmm
- make_truepheno
- make_phenotype
- make_gwas
- make_gwaslogp
- make_genofreq

Most of these methods accept three arguments:
- a frame (except 'make_empty_frame'),
- a config (will be discussed below), and
- custom parameters

The "config" has misc supportive information about frames, needed in order to implement this functionality.
Check find_config.m for an example. You can either edit find_config.m, or directly put parts of it into your script.
You have to adapt find_config with your custom paths local to your machine or cluster.


Pre-generated data
------------------

This repository contains code that generates realistic synthetic phenotypes.
Input data required to run this code can be accessed via `/projects/NS9114K/SYNGP` at NORSTORE (`ssh login.norstore.uio.no`).

The `/projects/NS9114K/SYNGP/` folder is organized as follows:

* `2558411_ref.bim` is a reference file containing 2558411 SNPs (also known as 2.5 million SNPs template).
* Folder `EUR_10K_2M_merged/` contains plink genotypes for 10K individuals and 2518104 SNPs (subset of the 2.5M template)
  This genotypes were generated with HAPGEN2 tool using EUR population from 1000 genome phase 3 as reference
* Folder `EUR_100K_80M_frames` contains pre-generated synthetic data calculated with `EUR_100K_80M` frame. Example:
  * `truepheno_<ID>_pi1e-03.mat` --- true phenotypes calculated from the model where `0.1%` of SNPs were causal (non-zero effect size), and the remaining SNPs had  effect size set to zero. `ID` is a unique identifier of the run, used to differentiate multiple runs with the same parameter.
  * `genofreq.mat` --- matrix of genotype frequencies for `EUR_100K_80M` frame.
* Folder `EUR_100K_1M_frames` contains pre-generated synthetic data for `EUR_100K_1M` frame.
  * `gwaslogp_<ID>_pi=1e-03_h2=0.90.mat` contains gwas results of the phenotypes. True phenotypes for this run were calculated using `EUR_100K_80M` frame, but GWAS results were conducted only using 1M snps.

* `annotations_RData/` --- annotation cathegory in `R` format. This is not yet integrated with the rest of the code.
