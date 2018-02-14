if ~exist('config', 'var'), find_config; end;         % => config

frame = make_empty_frame (config, 'EUR_10K_1188K');
frame = make_truebeta_gmm(frame, config, 'pivec', 0.0001, 'sig1vec', 1.0, 'sig2vec', 1.0, 'rhovec', 0.5);
frame = make_truepheno   (frame, config);
frame = make_phenotype   (frame, config, 'h2', 0.8);
frame = make_gwas        (frame, config);
frame = make_gwaslogp    (frame, config);

frame2 = reframe(frame, config, 'EUR_10K_1190K');
