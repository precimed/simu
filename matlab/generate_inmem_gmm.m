index     %index = 1;     % to be defined by caller

config = struct();

% A subset of 80M and 1M template (w_hm3)
config.frames.EUR_10K_1188K_merged.subj = 10000;
config.frames.EUR_10K_1188K_merged.snps = 1188973;
config.frames.EUR_10K_1188K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 10000, 'H:\NORSTORE\SYNGP\EUR_10K_1M_merged\all');
config.frames.EUR_10K_1188K_merged.bim = PlinkRead_bim('H:\NORSTORE\SYNGP\EUR_10K_1M_merged\all');

% Matlab 1M reference template (w_hm3)
config.frames.EUR_10K_1190K_ref.subj = 10000;
config.frames.EUR_10K_1190K_ref.snps = 1190321;
config.frames.EUR_10K_1190K_ref.bim  = PlinkRead_bim('H:\NORSTORE\SYNGP\1m_ref', true, '%s %s %f %d %s %s %s %s');

i = 1;
if index == i, pivec = []; sig1vec = []; sig2vec = []; rhovec = []; end; i = i + 1;  % no causal SNPs
if index == i, pivec = [0.0001]; sig1vec = [5]; sig2vec = [0]; rhovec = [0]; end; i = i + 1;  % causal only 1st phenotype
if index == i, pivec = [0.0001]; sig1vec = [5]; sig2vec = [5]; rhovec = [0]; end; i = i + 1;  % pleiotropy but no geneetic correlation 
if index == i, pivec = [0.0001]; sig1vec = [5]; sig2vec = [5]; rhovec = [0.98]; end; i = i + 1;  % pleiotropy
if index == i, pivec = [0.0001 0.0001]; sig1vec = [5 0.0]; sig2vec = [0.0 5]; rhovec = [0 0]; end; i = i + 1;  % indep
if index == i, pivec = [0.0001 0.0001 0.0001]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0.99]; end; i = i + 1;   % full
if index == i, pivec = [0.0001 0.0001]; sig1vec = [5 5]; sig2vec = [0.0 5]; rhovec = [0 0.99]; end; i = i + 1;   % first + pleio
if index == i, pivec = [0.0001 0.0001]; sig1vec = [5 5]; sig2vec = [5 5]; rhovec = [-0.99 0.99]; end; i = i + 1;   % pleio (opposite)

if index == i, pivec = [0.001]; sig1vec = [5]; sig2vec = [0]; rhovec = [0]; end; i = i + 1;  % causal only 1st phenotype
if index == i, pivec = [0.001]; sig1vec = [5]; sig2vec = [5]; rhovec = [0]; end; i = i + 1;  % pleiotropy but no geneetic correlation 
if index == i, pivec = [0.001]; sig1vec = [5]; sig2vec = [5]; rhovec = [0.98]; end; i = i + 1;  % pleiotropy
if index == i, pivec = [0.001 0.001]; sig1vec = [5 0.0]; sig2vec = [0.0 5]; rhovec = [0 0]; end; i = i + 1;  % indep
if index == i, pivec = [0.001 0.001 0.001]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0.99]; end; i = i + 1;   % full
if index == i, pivec = [0.001 0.001]; sig1vec = [5 5]; sig2vec = [0.0 5]; rhovec = [0 0.99]; end; i = i + 1;   % first + pleio
if index == i, pivec = [0.001 0.001]; sig1vec = [5 5]; sig2vec = [5 5]; rhovec = [-0.99 0.99]; end; i = i + 1;   % pleio (opposite)

if index == i, pivec = [0.01]; sig1vec = [5]; sig2vec = [0]; rhovec = [0]; end; i = i + 1;  % causal only 1st phenotype
if index == i, pivec = [0.01]; sig1vec = [5]; sig2vec = [5]; rhovec = [0]; end; i = i + 1;  % pleiotropy but no geneetic correlation
if index == i, pivec = [0.01]; sig1vec = [5]; sig2vec = [5]; rhovec = [0.98]; end; i = i + 1;  % pleiotropy
if index == i, pivec = [0.01 0.01]; sig1vec = [5 0.0]; sig2vec = [0.0 5]; rhovec = [0 0]; end; i = i + 1;  % indep
if index == i, pivec = [0.01 0.01 0.01]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0.99]; end; i = i + 1;   % full
if index == i, pivec = [0.01 0.01]; sig1vec = [5 5]; sig2vec = [0.0 5]; rhovec = [0 0.99]; end; i = i + 1;   % first + pleio
if index == i, pivec = [0.01 0.01]; sig1vec = [5 5]; sig2vec = [5 5]; rhovec = [-0.99 0.99]; end; i = i + 1;   % pleio (opposite)

i = 101;
if index == i, pivec = [0.0001 0.0001 0.0001]; sig1vec = [5 0.0 3]; sig2vec = [0.0 5 3]; rhovec = [0 0 0]; end; i = i + 1;  % indep (big) + pleio_no_correlation (small)
if index == i, pivec = [0.0001 0.0001 0.0001]; sig1vec = [3 0.0 5]; sig2vec = [0.0 3 5]; rhovec = [0 0 0]; end; i = i + 1;  % indep (small) + pleio_no_correlation (big)
if index == i, pivec = [0.0001 0.0001 0.0001 0.0001]; sig1vec = [5 0.0 3 3]; sig2vec = [0.0 5 3 3]; rhovec = [0 0 -0.99 0.99]; end; i = i + 1;   % indep (big) + pleio opposite (small) 
if index == i, pivec = [0.0001 0.0001 0.0001 0.0001]; sig1vec = [3 0.0 5 5]; sig2vec = [0.0 3 5 5]; rhovec = [0 0 -0.99 0.99]; end; i = i + 1;   % indep(small) pleio (opposite), big

if index == i, pivec = [0.001 0.001 0.001]; sig1vec = [5 0.0 3]; sig2vec = [0.0 5 3]; rhovec = [0 0 0]; end; i = i + 1;  % indep (big) + pleio_no_correlation (small)
if index == i, pivec = [0.001 0.001 0.001]; sig1vec = [3 0.0 5]; sig2vec = [0.0 3 5]; rhovec = [0 0 0]; end; i = i + 1;  % indep (small) + pleio_no_correlation (big)
if index == i, pivec = [0.001 0.001 0.001 0.001]; sig1vec = [5 0.0 3 3]; sig2vec = [0.0 5 3 3]; rhovec = [0 0 -0.99 0.99]; end; i = i + 1;   % indep (big) + pleio opposite (small) 
if index == i, pivec = [0.001 0.001 0.001 0.001]; sig1vec = [3 0.0 5 5]; sig2vec = [0.0 3 5 5]; rhovec = [0 0 -0.99 0.99]; end; i = i + 1;   % indep(small) pleio (opposite), big

if index == i, pivec = [0.01 0.01 0.01]; sig1vec = [5 0.0 3]; sig2vec = [0.0 5 3]; rhovec = [0 0 0]; end; i = i + 1;  % indep (big) + pleio_no_correlation (small)
if index == i, pivec = [0.01 0.01 0.01]; sig1vec = [3 0.0 5]; sig2vec = [0.0 3 5]; rhovec = [0 0 0]; end; i = i + 1;  % indep (small) + pleio_no_correlation (big)
if index == i, pivec = [0.01 0.01 0.01 0.01]; sig1vec = [5 0.0 3 3]; sig2vec = [0.0 5 3 3]; rhovec = [0 0 -0.99 0.99]; end; i = i + 1;   % indep (big) + pleio opposite (small) 
if index == i, pivec = [0.01 0.01 0.01 0.01]; sig1vec = [3 0.0 5 5]; sig2vec = [0.0 3 5 5]; rhovec = [0 0 -0.99 0.99]; end; i = i + 1;   % indep(small) pleio (opposite), big

frame = make_empty_frame (config, 'EUR_10K_1188K_merged');
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', sig1vec, 'sig2vec', sig2vec, 'rhovec', rhovec);
%plot(frame.truebeta(:, 1), frame.truebeta(:, 2), '.')

frame = make_truepheno   (frame, config, 'snpstep', 100);
frame = make_phenotype   (frame, config, 'h2', 0.9);

nsubj = frame.subj;
s500 = floor(nsubj/2); s501 = s500+1; s1000 = nsubj;
s333 = floor(nsubj/3); s334 = s333+1; s667 = 2*s333+1;

fprintf('n = no overlap\n');
frameN = make_gwas    (frame,  config, 'snpstep', 100, 'subj_idx', {1:s500, s501:s1000});
frameN = make_gwaslogp(frameN, config);
frameN = reframe      (frameN, config, 'EUR_10K_1190K_ref');

fprintf('p = partial overlap\n');
frameP = make_gwas    (frame,  config, 'snpstep', 100, 'subj_idx', {1:s667, s334:s1000});
frameP = make_gwaslogp(frameP, config);
frameP = reframe      (frameP, config, 'EUR_10K_1190K_ref');

fprintf('f = full overlap\n');
frameF = make_gwas    (frame,  config, 'snpstep', 100, 'subj_idx', {1:s1000, 1:s1000});
frameF = make_gwaslogp(frameF, config);
frameF = reframe      (frameF, config, 'EUR_10K_1190K_ref');

save(sprintf('result_%i.mat', index), 'frameF', 'frameN', 'frameP');

%plot(betavec(:, 1), betavec(:, 2), '.')
%plot(zvec1n, zvec2n, '.')
%plot(zvec1p, zvec2p, '.')
%plot(zvec1f, zvec2f, '.')




