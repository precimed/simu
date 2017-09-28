% Data folder with synthetic genotypes (can be downloaded from dropbox)
% https://www.dropbox.com/sh/b2it1oax6l0nc9b/AACT6ADa6bcr2gUV8tRwuvdka?dl=0
datafolder = 'H:\Dropbox\analysis\2017_09_September_28_Norbis';

% System code - bootstrap simulations...
if ~exist('config', 'var')
    fprintf('Setup genotypes... ');
    config.frames.EUR_10K_1188K_merged.subj = 10000;
    config.frames.EUR_10K_1188K_merged.snps = 1188973;
    config.frames.EUR_10K_1188K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 10000, fullfile(datafolder, 'EUR_10K_1M_merged/all'));
    config.frames.EUR_10K_1188K_merged.bim = PlinkRead_bim(fullfile(datafolder, 'EUR_10K_1M_merged/all'));
end
frame = make_empty_frame (config, 'EUR_10K_1188K_merged');

% Simulation - Step1. Generate additive effect sizes for 2 traits using
% 3-component mixture: first components has 5e-4 causal variants,
% affecting both traits, with 0.9 correlation of effect sizes.
% Second component is the same but correlation is set to -0.9.
frame = make_truebeta_gmm(frame, config, 'pivec', [5e-4 5e-4], 'sig1vec', [1 1], 'sig2vec', [1 1], 'rhovec', [0.9 -0.9]);
% now frame.truebeta gives a matrix of size [nsnp, 2] with true (causal) effect
% sizes in two simulated traits for each SNP

% Calculate genetic component of the phenotypes
frame = make_truepheno   (frame, config, 'snpstep', 100);

% Calculate actual phenotypes, setting heritability to h2=0.5
frame = make_phenotype   (frame, config, 'h2', 0.5);
% now frame.phenotype gives a matrix of size [nsubj, 2]
% with quantitative phenotype for each of the two traits

% Perform GWAS using non-overlapping subjects (5K subjects for each trait)
frame = make_gwas        (frame, config, 'subj_idx', {1:5000, 5001:10000});
% now frame.gwaspval and frame.gwasbeta contains estimated loci effect
% sizes and association p-value for SNP.

save('results', 'frame')
% The order of SNPs in the frame is the same as in the input .bim file.
