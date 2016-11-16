# SynGP
Synthetic genotypes and phenotypes

This repository contains code that generates realistic synthetic phenotypes.
Input data required to run this code can be accessed via `/projects/NS9114K/SYNGP` at NORSTORE (`ssh login.norstore.uio.no`).

The `/projects/NS9114K/SYNGP/` folder is organized as follows:

* `2558411_ref.bim` is a reference file containing 2558411 SNPs (also known as 2.5 million SNPs template).
* Folder `EUR_10K_2M_merged/` contains plink genotypes for 10K individuals and 2518104 SNPs (subset of the 2.5M template)
  This genotypes were generated with HAPGEN2 tool using EUR population from 1000 genome phase 3 as reference
* Folder `EUR_10K_2M_pheno/` contains pre-generated synthetic phenotypes. The naming scheme is like this

   * `pi<N>_h<H>_iter<K>.mat`, where `N = 0..6` indicated that a fraction of `10^-N` SNPs are causals; `H = 0, 20, 50, 80` is heritability of the phenotype (in percentage); `K` is an iteration (each data is randomly generated several times)
   * `pi<N>_h<H>_iter<K>_enrich10at4.mat` --- same as before, but with differential enrichment; this means that causal SNPs are randomly distributed across certain annotation cathegory

  Each of these files will containing
    * `betavec` --- true effect sizes (assigned according to the model described above)
    * `logpvec` --- `-log10(pval)` where `pval` simulate GWAS results for given phenotype`
    * `zvec` --- signed `z` score that corresponds to `logpvec`
    * `p1`, `h2` --- true  proportion of causal SNPs and heritability

* `annotations_RData/` --- annotation cathegory in `R` format. This is not yet integrated with the rest of the code.
