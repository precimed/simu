TBD:

* explain model and refer to http://cnsgenomics.com/software/gcta/#GWASSimulation
* explain file formats, both input and output
* explain difference between ``--gcta-sigma`` and ``--norm-effect``
* explain causal regions feature, and why they must not overlap (this isn't our view of biology - but makes simulations simpler)
* explain how to compile or download pre-generated binary

```
SIMU v0.9.0 - library for simulation of GWAS summary statistics:
  -h [ --help ]             produce this help message
  --bfile arg               prefix for plink .bed/.bim/.fam file
  --bfile-chr arg           same as --bfile, but will automatically concatenate
                            .bed/.bim/.fam files split across 22 chromosomes.
                            If the filename prefix contains the symbol @, SIMU
                            will replace the @ symbol with chromosome numbers.
                            Otherwise, SIMU will append chromosome numbers to
                            the end of the filename prefix.
  --qt                      simulate quantitative trait
  --cc                      simulate case/control trait
  --num-traits arg (=1)     number of traits (either 1 or 2 traits are
                            supported)
  --k arg                   prevalence for case/control traits, by default 0.1;
                            one value per trait
  --ncas arg                number of cases, by default N*k, where N is sample
                            size; one value per trait
  --ncon arg                number of controls, by default N*(1-k), where N is
                            sample size; one value per trait
  --hsq arg                 heritability, by default 0.7; one value per trait
  --num-components arg (=1) number of components in the mixture
  --causal-pi arg           proportion of causal variants; by default 0.001;
                            one value per mixture component
  --causal-n arg            number of causal variants (alternative to
                            --causal-pi option); one value per mixture
                            component
  --causal-variants arg     file with a list of causal variants and,
                            optionally, their effect sizes; one file per
                            mixture component. This is an alternative option to
                            --causal-pi and --causal-n. See README.md file for
                            detailed description of file formats.
  --causal-regions arg      file with a list of non-overlapping regions to
                            distribute causal variants. Regions must be defined
                            as a list of variant names (e.g. RS numbers, one
                            value per line). Supplements --causal-pi or
                            --causal-n options; can not be used together with
                            --causal-variants option. One file per mixture
                            component
  --trait1-sigsq arg        variance of effect sizes for trait1 per causal
                            marker; by default 1.0; one value per mixture
                            component
  --trait2-sigsq arg        variance of effect sizes for trait2 per causal
                            marker; by default 1.0; one value per mixture
                            component
  --gcta-sigma              draw effect sizes with variance inversely
                            proportional to sqrt(2*p(1-p)), where p is allele
                            frequency.
  --norm-effect             report effect sizes w.r.t. normalized genotypes
                            (e.i. additively coded 0,1,2 genotypes devided by
                            sqrt(2*p(1-p)), where p is allele frequency).
                            Default behavior without --norm-effect is to report
                            effect size w.r.t. additively coded 0,1,2
                            genotypes.
  --rg arg                  coefficient of genetic correlation; by default 0.0;
                            one value per mixture component
  --seed arg                seed for random numbers generator (default is
                            time-dependent seed)
  --out arg (=simu)         prefix of the output files; will generate .pheno
                            file containing synthesized phenotypes; and
                            .*.causals files (one file per trait) containing
                            lists of causal variants and their effect sizes for
                            each component in the mixture. See README.md file
                            for detailed description of file formats.

Examples:

* Simulate a quantitative trait with 1% causal markers, with the heritability of 0.5:
  simu --bfile test --qt --causal-pi 0.01 --hsq 0.5

* Simulate a quantitative trait as above, using simulation model from GCTA GWAS Simulation:
  simu --bfile test --qt --causal-pi 0.01 --hsq 0.5 --gcta-sigma --norm-effect

* Simulate two quantitative traits with genetic correlation of 0.8 and the heritabilities of 0.2 and 0.6:
  simu --bfile test --qt --causal-pi 0.01 --num-traits 2 --hsq 0.2 0.6 --rg 0.8

* Simulate 500 cases and 500 controls with the heritability of liability of 0.5 and disease prevalence of 0.1:
  simu --bfile test --cc --k 0.1 --ncas 500 --ncon 500 --causal-pi 0.01 --hsq 0.5

* Simulate two quantitative traits with three genetic components:
      (1) causal variants specific to the first trait,
      (2) causal variants specific to the second trait,
      (3) shared causal variants with correlation of 0.8:
  simu --bfile test --qt --num-traits 2 --hsq 0.2 0.6  --num-components 3 \
       --causal-pi 0.01 0.01 0.01 --trait1-sigsq 1 0 1 --trait2-sigsq 0 1 1 --rg 0 0 0.8

* Simulate a quantitative trait with specific location of causal markers:
  simu --bfile test --qt --causal-variants causal.snplist --hsq 0.5

* Simulate a quantitative trait where n=10 causal markers are randomly distributed among specific region:
  simu --bfile test --qt --causal-n 10 --causal-regions snplist --hsq 0.5
```
