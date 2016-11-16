if strcmp(getenv('COMPUTERNAME'), 'DESKTOP-696SVB3'), simupop = 'H:\ftp\ip113.ucsd.edu\SimuPop'; end;
if ~exist('simupop', 'var'), simupop = '/space/syn03/1/data/oleksandr/SimuPop'; end;

if ~exist('annomat', 'var'), load('H:\ftp\ip113.ucsd.edu\old_25_SNPs_processed_summary\GWAS_Annot\annomat.mat'); end;

dataDir = fullfile(simupop, 'EUR_10K_2M_merged');
outputDir = fullfile(simupop, 'EUR_10K_2M_pheno');
refFile = fullfile(simupop, '2558411_ref');

% Align to 2558411 reference template; 2518104 SNPs overlap (98.42%)
bim_ref = PlinkRead_bim(refFile, true, '%s %s %f %d %s %s %s %s');
bim     = PlinkRead_bim(fullfile(dataDir, 'all'));
[in_ref, index_to_ref] = ismember(bim.snplist, bim_ref.snplist);
[in_index, ref_to_index] = ismember(bim_ref.snplist, bim.snplist);
if ~all(in_ref), error('Unknown SNPs detected in the data'); end;

nsubj    = 10000;
nsnp     = length(bim.snplist);
nsnp_ref = length(bim_ref.snplist);

G = @(snps)PlinkRead_binary2(nsubj,snps,fullfile(dataDir, 'all'));  % int8, nsubj * nsnp

pi1vec  = logspace(-6, 0, 7);       % [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.0];
h2vec   = [0.0 0.2 0.5 0.8];
Nfrac = 1.0;
snpstep = 10000;

for iter = 1:5
for pi1 = pi1vec
for h2 = h2vec

descr = sprintf('pi%i_h%i_iter%i_enrich10at4', -log10(pi1), 100*h2, iter);
fprintf('Running %s...\n ', descr);

nsnpCausal = max(floor(pi1*nsnp), 1);
enriched = in_index & annomat(:, 10) > 4.0;
if sum(enriched) < nsnpCausal, continue; end;
betavec_tmp = zeros(nsnp, 1);
betavec_tmp(ref_to_index(randsample(find(enriched), nsnpCausal))) = randn(nsnpCausal, 1);
[~, rhovec_tmp, pvec_tmp] = SynGP_inmem(G, nsubj, nsnp, [], h2, Nfrac, snpstep, betavec_tmp);
%[betavec_tmp, rhovec_tmp, pvec_tmp] = SynGP_inmem(G, nsubj, nsnp, pi1, h2, Nfrac, snpstep);

betavec = nan(nsnp_ref, 1); rhovec = nan(nsnp_ref, 1); pvec = nan(nsnp_ref, 1);
betavec(index_to_ref) = betavec_tmp;
rhovec(index_to_ref) = rhovec_tmp;
pvec(index_to_ref) = pvec_tmp;

logpvec = -log10(pvec);
logpvec(isinf(logpvec)) = -log10(realmin);
zvec = -norminv(10.^-logpvec/2).*sign(rhovec);

file = fullfile(outputDir, descr);
save(file, 'betavec', 'logpvec', 'zvec', 'pi1', 'Nfrac', 'h2', '-v7.3');
fprintf('Results are saved to %s.mat\n', file);

end
end
end
