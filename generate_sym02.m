pi_index        %pi_index = 1 .. 12;     % to be defined by caller
pe_index        %pe_index = 1 .. 4;

rng('shuffle')

if ~exist('config', 'var'), find_config; end;
fieldnames(config.frames)

pleiotropic_enrichment = { '_Independent', '_Correlation', '_SharedSnps', '_Opposite' };
overlap_options = { '_NoneOverlap', '_PartialOverlap', '_FullOverlap' };

hvec_opts =  [0.1, 0.3, 0.5, 0.7, 0.9];
pivec_opts = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1];

pi_index = mod(pi_index - 1, length(pivec_opts)) + 1;
if pi_index <= 0 || pi_index > length(pivec_opts), error('invalid pi_index'); end;
pival = pivec_opts(pi_index);

pe_index = mod(pe_index - 1, length(pleiotropic_enrichment)) + 1;
if pe_index <= 0 || pe_index > length(pleiotropic_enrichment), error('invalid pe_index'); end;

if pe_index == 1, pivec = [pival/2 pival/2]; sig1vec = [5 0.0]; sig2vec = [0.0 5]; rhovec = [0 0]; end;  % indep
if pe_index == 2, pivec = [pival/3 pival/3 pival/3]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0.99]; end;   % full
if pe_index == 3, pivec = [pival/3 pival/3 pival/3]; sig1vec = [5 0.0 5]; sig2vec = [0.0 5 5]; rhovec = [0 0 0]; end;   % indep (big) + pleio_no_correlation (small)
if pe_index == 4, pivec = [pival/4 pival/4 pival/4 pival/4]; sig1vec = [5 0.0 5 5]; sig2vec = [0.0 5 5 5]; rhovec = [0 0 -0.99 0.99]; end;   % indep (big) + pleio opposite (small) 

frame = make_empty_frame (config, 'EUR_100K_80M_chunks');
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', sig1vec, 'sig2vec', sig2vec, 'rhovec', rhovec);
frame = make_truepheno   (frame, config, 'snpstep', 100);

[~,filename,~] = fileparts(tempname);
save(sprintf('truepheno_%s_pi%.0e%s.mat', filename, pival, pleiotropic_enrichment{pe_index}), 'frame', '-v7.3');

nsubj = frame.subj;
s500 = floor(nsubj/2); s501 = s500+1; s1000 = nsubj;
s333 = floor(nsubj/3); s334 = s333+1; s667 = 2*s333+1;
no_overlap = {1:s500, s501:s1000}; % no overlap
partial_overlap = {1:s667, s334:s1000}; % partial overlap
full_overlap = {1:s1000, 1:s1000};   % full overlap
overlap_indices = {no_overlap, partial_overlap, full_overlap};

frame0 = frame;
for hindex = 1:length(hvec_opts)
    frame = reframe       (frame0, config, 'EUR_100K_1188K_merged');
    frame = make_phenotype(frame, config, 'h2', hvec_opts(hindex));
    for oindex = 1:length(overlap_options)
        frame = make_gwas     (frame, config, 'subj_idx', overlap_indices{oindex});
        frame = make_gwaslogp (frame, config);
        frame = reframe       (frame, config, 'EUR_100K_1190K_ref');
        frame.opts.now = datestr(now);
        save(sprintf('gwaslogp_%s_pi=%.0e_h2=%.2f%s%s.mat', filename, pival, hvec_opts(hindex), pleiotropic_enrichment{pe_index}, overlap_options{oindex}), 'frame', '-v7.3');
    end
end

