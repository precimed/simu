index 
%index = 1;     % to be defined by caller
%simupop = 'C:\Users\oleksanf\Dropbox\analysis\2016_09_September_22_Bivariate\SYNGP';
simupop = '/work/users/oleksanf';

reffile = fullfile(simupop, '1m_ref');
bfile = fullfile(simupop, 'EUR_10K_1M_merged', 'all');

nsubj    = 100000;
nsnp     = 1188973;
snpstep  = 100;

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

[betavec, mix] = make_gmm_beta(nsnp, pivec, sig1vec, sig2vec, rhovec);
%nsnp2 = floor(nsnp/2);betavec((nsnp2:nsnp)', 1) = 0;betavec((1:nsnp2)', 2) = 0;

%figure(index)
%plot(betavec(:, 1), betavec(:, 2), '.')

h2 = 0.9;

s500 = floor(nsubj/2); s501 = s500+1; s1000 = nsubj;
s333 = floor(nsubj/3); s334 = s333+1; s667 = 2*s333+1;

% n = no overlap
[zvec1n, nvec1n, betavec1n, rhovec1n, logpvec1n, mix1n] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec(:, 1), h2, 1:s500, snpstep, mix);
[zvec2n, nvec2n, betavec2n, rhovec2n, logpvec2n, mix2n] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec(:, 2), h2, s501:s1000, snpstep, mix);

% p = partial overlap
[zvec1p, nvec1p, betavec1p, rhovec1p, logpvec1p, mix1p] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec(:, 1), h2, 1:s667, snpstep, mix);
[zvec2p, nvec2p, betavec2p, rhovec2p, logpvec2p, mix2p] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec(:, 2), h2, s334:s1000, snpstep, mix);

% f = full overlap
[zvec1f, nvec1f, betavec1f, rhovec1f, logpvec1f, mix1f] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec(:, 1), h2, 1:s1000, snpstep, mix);
[zvec2f, nvec2f, betavec2f, rhovec2f, logpvec2f, mix2f] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec(:, 2), h2, 1:s1000, snpstep, mix);

save(sprintf('result_%i.mat', index));

%plot(betavec(:, 1), betavec(:, 2), '.')
%plot(zvec1n, zvec2n, '.')
%plot(zvec1p, zvec2p, '.')
%plot(zvec1f, zvec2f, '.')
