%pi_index        %pi_index = 1 .. 12;     % to be defined by caller
%pg_index        %pg_index = 1 .. 6;

rng('shuffle')

if ~exist('config', 'var'), find_config; end;
fieldnames(config.frames)

hvec_opts =  [0.9]; %[0.1, 0.3, 0.5, 0.7, 0.9];
pivec_opts = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1];
gencorr_opts = [NaN, 0, 0.5, 0.5, 0.75, 1.0];

pi_index = mod(pi_index - 1, length(pivec_opts)) + 1;
if pi_index <= 0 || pi_index > length(pivec_opts), error('invalid pi_index'); end;
pival = pivec_opts(pi_index);

pg_index = mod(pg_index - 1, length(gencorr_opts)) + 1;
if pg_index <= 0 || pg_index > length(gencorr_opts), error('invalid pg_index'); end;
gencorr = gencorr_opts(pg_index);

if isnan(gencorr), pivec = [pival, pival]; sig1vec = [1 0]; sig2vec = [0 1]; rhovec = [0 0];  % indep
else pivec = pival; sig1vec = [1]; sig2vec = [1]; rhovec = gencorr;                           % pleio
end

frame = make_empty_frame (config, 'EUR_100K_80M_chunks');
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', sig1vec, 'sig2vec', sig2vec, 'rhovec', rhovec);
frame = make_truepheno   (frame, config, 'snpstep', 100);

[~,filename,~] = fileparts(tempname);
save(sprintf('truepheno_%s_pi=%.0e_gencorr=%.2f.mat', filename, pival, gencorr), 'frame', '-v7.3');

nsubj = frame.subj;
s500 = floor(nsubj/2); s501 = s500+1; s1000 = nsubj;
s333 = floor(nsubj/3); s334 = s333+1; s667 = 2*s333+1;
no_overlap = {1:s500, s501:s1000}; % no overlap
partial_overlap = {1:s667, s334:s1000}; % partial overlap
full_overlap = {1:s1000, 1:s1000};   % full overlap
overlap_indices = {no_overlap, partial_overlap, full_overlap};

frame0 = frame;
for hindex = 1:length(hvec_opts)
    frame1 = reframe       (frame0, config, 'EUR_100K_8801K_merged');
    frame1 = make_phenotype(frame1, config, 'h2', hvec_opts(hindex));
    frame = make_gwas     (frame1, config, 'subj_idx', no_overlap);
    frame = make_gwaslogp (frame, config);
    frame = reframe       (frame, config, 'EUR_100K_9279K_ref');
    frame.opts.now = datestr(now);
    save(sprintf('gwaslogp_%s_pi=%.0e_gencorr=%.2f_h2=%.2f.mat', filename, pival, gencorr, hvec_opts(hindex)), 'frame', '-v7.3');
end

