% Frames could be of 3 distinct types:
% merged - a set of 3 files with all genotypes (all.bed, all.bim, all.fam)
% chunk  - chunked files, one per chromosome and position
% ref    - reference .bim file with columns CHR,SNP,GP,BP,A1,A2

config = struct();

% A subset of 80M  and 2.5M template
config.frames.EUR_10K_2518K_merged.snps = 2518104;
config.frames.EUR_10K_2518K_merged.subj = 10000;
config.frames.EUR_10K_2518K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 10000, 'H:\ftp\ip113.ucsd.edu\SimuPop\EUR_10K_2M_merged\all');

% Matlab 2.5M reference template
config.frames.EUR_10K_2558K_ref.subj = 10000;
config.frames.EUR_10K_2558K_ref.snps = 2558411;

% A subset of 80M and 1M template (w_hm3)
config.frames.EUR_10K_1188K_merged.subj = 10000;
config.frames.EUR_10K_1188K_merged.snps = 1188973;
config.frames.EUR_10K_1188K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 10000, 'H:\NORSTORE\SYNGP\EUR_10K_1M_merged\all');
config.frames.EUR_10K_1188K_merged.bim = PlinkRead_bim('H:\NORSTORE\SYNGP\EUR_10K_1M_merged\all');

% Matlab 1M reference template (w_hm3)
config.frames.EUR_10K_1190K_ref.subj = 10000;
config.frames.EUR_10K_1190K_ref.snps = 1190321;
config.frames.EUR_10K_1190K_ref.bim  = PlinkRead_bim('H:\NORSTORE\SYNGP\1m_ref', true, '%s %s %f %d %s %s %s %s');

% 80M chunked frame - 80378054 SNPs, 3920 non-empty chunks, 3931 chunks (11 weird empty chunks)
config.frames.EUR_100K_80M_chunks.subj = 100000;
config.frames.EUR_100K_80M_chunks.snps = 80378054;
config.frames.EUR_100K_80M_chunks.reader = ...
    @(subjvec, snpvec)chunked_frame_reader(subjvec, snpvec, ...
        100000, ...  % number of subjects
        'E:\EUR_100K_81M_gzip\misc\snp_detail.mat',...
        'E:\EUR_100K_81M_gzip\EUR_100K_81M_gzip');

% Set 'name' field for each frame
%fn = fieldnames(config.frames); for i=1:length(fn), config.frames.(fn{i}).name = fn{i}; end;