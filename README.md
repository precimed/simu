```
SIMU v0.9.0 - library for simulation of GWAS summary statistics:
  -h [ --help ]             produce help message
  --bfile arg               Prefix for Plink .bed/.bim/.fam file
  --bfile-chr arg           Same as --bfile, but will automatically concatenate
                            .bed/.bim/.fam files split across 22 chromosomes.
                            If the filename prefix contains the symbol @, SIMU
                            will replace the @ symbol with chromosome numbers.
                            Otherwise, SIMU will append chromosome numbers to
                            the end of the filename prefix.
  --qt                      simulate quantitative trait
  --cc                      simulate case/control trait
  --num-traits arg (=1)     Number of traits (either 1 or 2 traits are
                            supported)
  --k arg                   prevalence for case/control traits, by default 0.1;
                            one value per trait
  --ncas arg                number of cases, by default N*k; one value per
                            trait
  --ncon arg                number of controls, by default N*(1-k); one value
                            per trait
  --hsq arg                 heritability, by default 0.7; one value per trait
  --num-components arg (=1) Number of components in the mixture
  --causal-pi arg           proportion of causal variants; by default 0.001;
                            one value per mixture component
  --trait1-sigsq arg        variance of effect sizes for trait1 per causal
                            marker; by default 1.0; one value per mixture
                            component
  --trait2-sigsq arg        variance of effect sizes for trait2 per causal
                            marker; by default 1.0; one value per mixture
                            component
  --rg arg                  coefficient of genetic correlation; by default 0.0;
                            one value per mixture component
  --out arg (=simu)         Prefix of the output file; will generate .pheno
                            file (phenotypes) and .1.causals file (one per
                            trait, list MarkerName for all causal variants and
                            their effect sizes.
  --seed arg                Seed for random numbers generator (default is
                            time-dependent seed)
```
