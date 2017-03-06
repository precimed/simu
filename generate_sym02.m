index        %index = 1;     % to be defined by caller

rng('shuffle')

if ~exist('config', 'var'), find_config; end;
fieldnames(config.frames)

pleiotropic_enrichment = { '_Independent', '_Correlation', '_SharedSnps', '_Opposite' };

hvec_opts =  [0.1, 0.3, 0.5, 0.7, 0.9];
pivec_opts = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1];

index = mod(index - 1, length(pivec_opts)) + 1;
if index <= 0 || index > length(pivec_opts), error('invalid index'); end;
pival = pivec_opts(index);

for di = 1:4

if di == 1, pivec = [pival/2 pival/2]; sig1vec = [5 0.0]; sig2vec = [0.0 5]; rhovec = [0 0]; end;  % indep
if di == 2, pivec = [pival/3 pival/3 pival/3]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0.99]; end;   % full
if di == 3, pivec = [pival/3 pival/3 pival/3]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0]; end;   % indep (big) + pleio_no_correlation (small)
if di == 4, pivec = [pival/4 pival/4 pival/4 pival/4]; sig1vec = [5 0.0 5 5]; sig2vec = [0.0 5 5 5]; rhovec = [0 0 -0.99 0.99]; end;   % indep (big) + pleio opposite (small) 

frame = make_empty_frame (config, 'EUR_100K_80M_chunks');
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', sig1vec, 'sig2vec', sig2vec, 'rhovec', rhovec);
frame = make_truepheno   (frame, config, 'snpstep', 100);

[~,filename,~] = fileparts(tempname);
save(sprintf('truepheno_%s_pi%.0e%s.mat', filename, pival, pleiotropic_enrichment{di}), 'frame', '-v7.3');

frame0 = frame;
for hindex = 1:length(hvec_opts)
    frame = reframe       (frame0, config, 'EUR_100K_1188K_merged');
    frame = make_phenotype(frame, config, 'h2', hvec_opts(hindex));
    frame = make_gwas     (frame, config);
    frame = make_gwaslogp (frame, config);
    frame = reframe       (frame, config, 'EUR_100K_1190K_ref');
    frame.opts.now = datestr(now);
    save(sprintf('gwaslogp_%s_pi=%.0e_h2=%.2f%s.mat', filename, pival, hvec_opts(hindex), pleiotropic_enrichment{di}), 'frame', '-v7.3');
end

end
