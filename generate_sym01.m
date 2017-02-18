index        %index = 1;     % to be defined by caller

if ~exist('config', 'var'), find_config; end;
fieldnames(config.frames)

hvec_opts =  [0.1, 0.3, 0.5, 0.7, 0.9];
pivec_opts = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2];

index = mod(index - 1, length(hvec_opts)) + 1;
if index <= 0 || index > length(pivec_opts), error('invalid index'); end;
pivec = pivec_opts(index);

frame = make_empty_frame (config, 'EUR_100K_80M_chunks');
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', 1.0, 'sig2vec', 1.0, 'rhovec', 0.0);
        frame.truebeta = frame.truebeta(:, 1);
        frame.mixture  = frame.mixture(:, 1);
frame = make_truepheno   (frame, config, 'snpstep', 100);

[~,filename,~] = fileparts(tempname);
save(sprintf('truepheno_%s_pi%.0e.mat', filename, pivec), 'frame', '-v7.3');

%load(fullfile(norstore_path, 'SYNGP/EUR_100K_80M_pheno/polygenic_20170218_index5.mat'))
frame0 = frame;
for hindex = 1:length(hvec_opts)
    frame = reframe       (frame0, config, 'EUR_100K_1188K_merged');
    frame = make_phenotype(frame, config, 'h2', hvec_opts(hindex));
    frame = make_gwas     (frame, config);
    frame = make_gwaslogp (frame, config);
    frame = reframe       (frame, config, 'EUR_100K_1190K_ref');
    save(sprintf('gwaslogp_%s_pi=%.0e_h2=%.2f.mat', filename, pivec, hvec_opts(hindex)), 'frame', '-v7.3');
end

