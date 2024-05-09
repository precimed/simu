SIMU simulates a GWAS based on real genotype data.
This tool is a extension of ``gcta --simu-qt`` and ``gcta --simu-cc`` functionality,
described [here](http://cnsgenomics.com/software/gcta/#GWASSimulation).
Additional features include
* simulation of two traits (``--num-traits 2``), with given genetic correlation (``--rg``)
* allow input genotypes to be split by chromosome (``--bfile-chr``)
* efficient memory usage (most simulations require less than 2 GB of RAM; only ``.bim`` and ``.fam`` files are loaded into memory - ``.bed`` files are read in small chunks as needed)

``SIMU`` implements the same model as GCTA,
e.i. simple additive genetic model ``y=Gx+e``,
where `y` is output phenotype, ``G`` is genotype matrix, ``x`` is vector of effect sizes, ``e`` is environmental noise.

NB! Both ``SIMU`` and ``GCTA`` generate ``x`` from normal distribution.
There is one important difference here.
In ``GCTA``, effect sizes are in the units of normalized genotypes, e.i. 0,1,2 genotypes divided by ``sqrt(2*p(1-p))``, where p is allele frequency.
In ``SIMU``, by default, effect sizes are in the units of additively coded 0,1,2 genotypes.
Thus, the default behavior of ``SIMU`` is different from GCTA.
You may force ``SIMU`` to use ``GCTA`` model by specifying ``--gcta-sigma`` flag.
In this case you may also want to use ``--norm-effect`` flag to report effect sizes in the units of normalized genotypes.

Getting started
---------------

``SIMU`` is written in `C++` and is available for Linux / Unix and MacOS, but not for Windows.
You may download pre-compiled SIMU binary from the [Releases page](https://github.com/precimed/simu/releases) of this repository.
To test type ``simu --help``, which should produce a list of available options.
You may also build ``SIMU`` on your machine (see instructions further below).

To run ``SIMU`` you need raw genotypes.
You may download demo data from ``mostest_demo.tar.gz`` file from [here](https://1drv.ms/u/s!Ai1YZmdFa9ati5J246ZMEFWaUdZeFQ?e=FckenT), or a larger demo data from [here](https://1drv.ms/u/s!Ai1YZmdFa9atjIYmpZXUa_-XG2ur9Q?e=3oNkoK).

We've also prepared a set of synthetic genotypes for 100K ``individuals`` and ca. 11M markers.
It was produced by running [hapgen2](http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) software using EUR population from 1000 Genomes Phase 3 data as a reference (503 individuals). This set heavy and can be shared upon request.

Getting help
------------

If things didn't work, or if you have any suggestions, please submit [a new issue](https://github.com/precimed/simu/issues/new).
You pull requests are also very welcome!

File formats
------------

``simu`` reads genotypes from plink `.bed`/`.bim`/`.fam` files.

Output ``.pheno`` file contains tab-separated table, with header,
followed by one row per individual in the same order as the input `.fam` files.
The table contains three or four columnes.
First two columns are family and individual indeitifier from `.fam` file,
the remaining one or two columns contain simulated phenotypes.
For case/control traits,
``1`` means unaffected (control), 
``2`` means affected (case),
``-9`` means missing phenotype.
Example of ``.pheno`` file for quantitative trait (`--qt`):
```
FID     IID     trait1
id1_0   id2_0   -0.841613
id1_1   id2_1   -0.0747612
id1_2   id2_2   -2.04002
```

Output ``.causals`` file contains tab-separated file with effect sizes of causal variants.
It contains only markers with non-zero effect size, e.i. subset of variants from input ``.bim`` file.
First 5 columns contain information from ``.bim`` file:
``SNP`` (marker name or RS number),
``CHR`` (chromosome position),
``BP`` (base-pair positino),
``A1``, ``A2`` (reference and other allele).
The next column, ``FRQ``, contains frequency of ``A1`` allele.
Subsequent columns contain effect size w.r.t. ``A1``, one column for each of the ``--num-component`` components.
Depending on ``--norm-effect`` flag, effect sizes could be either in
units of additively coded 0,1,2 genotypes (default behavior), or
in units of normalized genotypes, e.i. 0,1,2 genotypes divided by ``sqrt(2*p(1-p))``, where p is allele frequency.
Example of ``.causals`` file:
```
SNP     CHR     POS     A1      A2      FRQ     BETA_c1
rs72651487      1       19706089        A       G       0.285715        -0.0335917
rs28631635      1       25728930        G       T       0.432055        0.0603244
rs67388349      1       28722149        C       CT      0.37678 -0.0527829
```

An optional argument ``--causal-variants`` allows to specify a file that explicitly lists causal variants and,
optionally, their effect sizes. The file must be whitespace-delimited table without header, containing one row per causal variant.
First column must contain marker name (to be matched with SNP column from `.bim` file). 
The remaining columns may contain effect sizes. For two-trait simulations there must be two effect size columns, e.i. effect size in each trait. Depending on ``--norm-effect`` flag, input effect sizes will be treated either as in
units of additively coded 0,1,2 genotypes (default behavior), or
in units of normalized genotypes, e.i. 0,1,2 genotypes divided by ``sqrt(2*p(1-p))``, where p is allele frequency.
An example of file for ``--causal-variants``:
```
rs72651487 -0.0335917
rs28631635 0.0603244
rs67388349 -0.0527829
```
The format for an optional argument ``--causal-regions`` is the same as ``--causal-variants``, except it must have only the list of marker names, but not their effect sizes.

Command line options
--------------------
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
  --trait1-s-pow arg        draw effect sizes on the first trait with variance
                            proportional to (2*p(1-p))^S, where parameter S is
                            defined by trait1-s-pow, and p is allele frequency;
                            by default 0.0; one value per mixture component
  --trait2-s-pow arg        draw effect sizes on the second trait with variance
                            proportional to (2*p(1-p))^S, where parameter S is
                            defined by trait2-s-pow, where p is allele
                            frequency; by default 0.0; one value per mixture
                            component
  --gcta-sigma              draw effect sizes with variance inversely
                            proportional to 2*p(1-p), where p is allele
                            frequency; this corresponds to --trait1-s-pow and
                            --trait2-s-pow set to -1.0
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
  --trait2-snp-offset arg   shifts causal variant in trait2 by a given
                            number of positions (for example,
                            '--trait2-snp-offset 1' indicates that causal
                            SNPs will correspond to adjacent rows in the
                            input bim file)
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

Compile SIMU on your machine
----------------------------

Check out [Releases page](https://github.com/precimed/simu/releases) for pre-compiled binaries.

* Step 1. Install prerequisites. On Ubuntu: ``sudo apt-get install git build-essential libboost-all-dev cmake``.
* Step 2. Clone this repository: ``git clone --recurse-submodules https://github.com/precimed/simu.git``.
  Note that this repository includes submodules. If you use older version of git you may want to review [this](https://stackoverflow.com/questions/3796927/how-to-git-clone-including-submodules).
* Step 3. Install to default location:
```
mkdir build && cd build
cmake ..
make && sudo make install
```
This will create place ``simu`` executable to ``/usr/local/bin`` folder.
* Step 3'. Install to a custom location:
```
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=~
make && make install
```
This will create place ``simu`` executable to your ``$HOME/bin`` folder.
* Step 4. Enjoy, ``simu`` is ready. Type ``simu --help`` to list available options.
