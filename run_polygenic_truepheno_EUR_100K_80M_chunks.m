index        %index = 1;     % to be defined by caller
framename    %framename = 'EUR_100K_80M_chunks';  % to be defined by caller
config = struct();

t = 'H:\NORSTORE';                                if exist(t, 'dir'), norstore_path = t; end;
t = '/work/users/oleksanf/NORSTORE';              if exist(t, 'dir'), norstore_path = t; end;
t = 'E:\EUR_100K_80M';                            if exist(t, 'dir'), path80m = t; end;
t = '/work/users/oleksanf/MMIL/cfan_SimuPop/EUR'; if exist(t, 'dir'), path80m = t; end;
snp_detail_path = fullfile(norstore_path, 'SYNGP/misc/snp_detail.mat');

% 80M chunked frame - 80378054 SNPs, 3920 non-empty chunks, 3931 chunks (11 weird empty chunks)
config.frames.EUR_100K_80M_chunks.subj = 100000;
config.frames.EUR_100K_80M_chunks.snps = 80378054;
config.frames.EUR_100K_80M_chunks.reader = @(subjvec, snpvec)chunked_frame_reader(subjvec, snpvec, 100000, snp_detail_path, path80m);

% 80M chunked frame, 500 subjects
config.frames.EUR_500_80M_chunks.subj = 500;
config.frames.EUR_500_80M_chunks.snps = 80378054;
config.frames.EUR_500_80M_chunks.reader = @(subjvec, snpvec)chunked_frame_reader(subjvec, snpvec, 500, snp_detail_path, fullfile(norstore_path, 'SYNGP/EUR_500_80M'));

%framename = 'EUR_500_80M_chunks';
%framename = 'EUR_100K_80M_chunks';

pivec_opts = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2];
if index <= 0 || index > length(pivec_opts), error('invalid index'); end;
pivec = pivec_opts(index);

frame = make_empty_frame (config, framename);
%frame = make_genofreq   (frame, config, 'snpstep', 10000); save(sprintf('genofreq_%s.mat', framename), 'frame', '-v7.3');
frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', 1.0, 'sig2vec', 1.0, 'rhovec', 0.0);
frame = make_truepheno   (frame, config, 'snpstep', 100);

[~,filename,~] = fileparts(tempname);
save(sprintf('result_polygenic_%i_%s_%s.mat', index, framename, filename), 'frame', '-v7.3');


