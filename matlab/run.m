% Simulate a pair of synthetic quantitative traits
% Supports the following features:
% pi        - polygenicity of both traits; default: 1e-3
% h2        - heritability of both traits; default: 0.5
% pi12      - polygenic overlap between traits; default: pi
%             limitation: must not exceed pi.
% rho12     - genetic correlation within polygenic overlap; default: 0
% n12       - number of overlaping subjects across GWASs
%             limitation: must not exceed nsubj (frame parameter)
% repeat    - how many traits to make with the same parameters
%             each time different SNPs will be causal, and different noise
%             applied to the phenotype
%
% Options that control which set of SNPs to use:
% (you may need up update paths in find_config.m).
%
% frame_pheno - frame to use when generate phenotype, default: EUR_100K_80M_chunks
% frame_gwas  - frame to use when perform gwas, default: EUR_100K_1188K_merged
% frame_final - frame to use when save the result, default: EUR_100K_1190K_ref
%
% Features available in massim but not supported by this script:
% differential_enrichment

rng('shuffle')

if ~exist('pi', 'var'), pi = 1e-3; end
if ~exist('h2', 'var'), h2 = 0.5; end
if ~exist('pi12', 'var'), pi12 = pi; end
if ~exist('rho12', 'var'), rho12 = 0; end
if ~exist('n12', 'var'), n12 = 0; end
if ~exist('repeat', 'var'), repeat = 1; end
if ~exist('frame_pheno', 'var'), frame_pheno = 'EUR_100K_80M_chunks'; end
if ~exist('frame_gwas', 'var'), frame_gwas = 'EUR_100K_1188K_merged'; end
if ~exist('frame_final', 'var'), frame_final = 'EUR_100K_1190K_ref'; end

if pi > 1 || pi < 0 || isnan(pi), error('pi out of range'); end;
if pi12 > 1 || pi12 < 0 || isnan(pi12) || pi12 > pi, error('pi12 out of range'); end;
if 2*pi - pi12 > 1, error('parameters out of range: pi + pi - pi12 > 1'); end;
pi1 = pi - pi12;  % now pi1 represents the weight of the component specific to each trait (excluding overlap)

if h2 < 0 || h2 > 1 || isnan(h2), error('h2 out of range'); end;
if rho12 < -1 || rho12 > 1 || isnan(rho12), error('rho12 out of range'); end;

if ~exist('config', 'var'), find_config; end;
fieldnames(config.frames)
disp(config.frames.(frame_pheno))
disp(config.frames.(frame_gwas))
disp(config.frames.(frame_final))

if n12 < 0 || n12 > config.frames.(frame_pheno).subj, error('n12 out of range'); end;

pivec = [];sig1vec = [];sig2vec = [];rhovec = [];
if pi1>0, pivec = [pivec pi1 pi1]; sig1vec = [sig1vec 1 0]; sig2vec = [sig2vec 0 1]; rhovec = [rhovec 0 0]; end;
if pi12>0, pivec = [pivec pi12]; sig1vec = [sig1vec 1]; sig2vec = [sig2vec 1]; rhovec = [rhovec rho12]; end;

for irepeat = 1:repeat
frame = make_empty_frame (config, frame_pheno);
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', sig1vec, 'sig2vec', sig2vec, 'rhovec', rhovec);
frame = make_truepheno   (frame, config, 'snpstep', 100);

nsubj = frame.subj;
subj_idx = {1:fix((nsubj+n12)/2), 1+fix((nsubj-n12)/2):nsubj};

frame = reframe       (frame, config, frame_gwas);
frame = make_phenotype(frame, config, 'h2', h2);
frame = make_gwas     (frame, config, 'subj_idx', subj_idx);
frame = make_gwaslogp (frame, config);
frame = reframe       (frame, config, frame_final);
frame.opts.now = datestr(now);

[~,filename,~] = fileparts(tempname);
save(sprintf('gwaslogp_h2=%.2f_pi=%.0e_pi12=%.0e_rho12=%.2f_n12=%i_%s.mat', h2, pi, pi12, rho12, fix(n12), filename), 'frame', '-v7.3');
end