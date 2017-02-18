% Frames could be of 3 distinct types:
% merged - a set of 3 files with all genotypes (all.bed, all.bim, all.fam)
% chunk  - chunked files, one per chromosome and position
% ref    - reference .bim file with columns CHR,SNP,GP,BP,A1,A2

config = struct();

t = 'H:\NORSTORE';                       if exist(t, 'dir'), norstore_path = t; end;
t = '/work/users/oleksanf/NORSTORE';     if exist(t, 'dir'), norstore_path = t; end;
t = 'E:\EUR_100K_80M';                   if exist(t, 'dir'), EUR_100K_80M_path = t; end;
t = '/work/users/oleksanf/EUR_100K_80M'; if exist(t, 'dir'), EUR_100K_80M_path = t; end;
snp_detail_path = fullfile(norstore_path, 'SYNGP/misc/snp_detail.mat');

% 80M chunked frame - 80378054 SNPs, 3920 non-empty chunks, 3931 chunks (11 weird empty chunks)
config.frames.EUR_100K_80M_chunks.subj = 100000;
config.frames.EUR_100K_80M_chunks.snps = 80378054;
config.frames.EUR_100K_80M_chunks.reader = @(subjvec, snpvec)chunked_frame_reader(subjvec, snpvec, 100000, snp_detail_path, EUR_100K_80M_path);

% 80M chunked frame, 10K subjects (reader is not available yet)
%config.frames.EUR_10K_80M_chunks.subj = 10000;
%config.frames.EUR_10K_80M_chunks.snps = 80378054;

% 80M chunked frame, 500 subjects
%config.frames.EUR_500_80M_chunks.subj = 500;
%config.frames.EUR_500_80M_chunks.snps = 80378054;
%config.frames.EUR_500_80M_chunks.reader = @(subjvec, snpvec)chunked_frame_reader(subjvec, snpvec, 500, snp_detail_path, fullfile(norstore_path, 'SYNGP/EUR_500_80M'));

% A subset of 80M  and 2.5M template
%config.frames.EUR_10K_2518K_merged.snps = 2518104;
%config.frames.EUR_10K_2518K_merged.subj = 10000;
%config.frames.EUR_10K_2518K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 10000, fullfile(norstore_path, 'SYNGP/EUR_10K_2M_merged/all'));
%config.frames.EUR_10K_2518K_merged.bim = PlinkRead_bim(fullfile(norstore_path, 'SYNGP/EUR_10K_2M_merged/all'));

% Matlab 2.5M reference template
%config.frames.EUR_10K_2558K_ref.subj = 10000;
%config.frames.EUR_10K_2558K_ref.snps = 2558411;
%config.frames.EUR_10K_2558K_ref.bim  = PlinkRead_bim(fullfile(norstore_path, 'SYNGP/2558411_ref'), true, '%s %s %f %d %s %s %s %s');

% A subset of 80M and 1M template (w_hm3) -- 10K subjects
%config.frames.EUR_10K_1188K_merged.subj = 10000;
%config.frames.EUR_10K_1188K_merged.snps = 1188973;
%config.frames.EUR_10K_1188K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 10000, fullfile(norstore_path, 'SYNGP/EUR_10K_1M_merged/all'));
%config.frames.EUR_10K_1188K_merged.bim = PlinkRead_bim(fullfile(norstore_path, 'SYNGP/EUR_10K_1M_merged/all'));

% Matlab 1M reference template (w_hm3) -- 10K subjects
%config.frames.EUR_10K_1190K_ref.subj = 10000;
%config.frames.EUR_10K_1190K_ref.snps = 1190321;
%config.frames.EUR_10K_1190K_ref.bim  = PlinkRead_bim(fullfile(norstore_path, 'SYNGP/1m_ref'), true, '%s %s %f %d %s %s %s %s');

% A subset of 80M and 1M template (w_hm3) -- 100K subjects
config.frames.EUR_100K_1188K_merged.subj = 100000;
config.frames.EUR_100K_1188K_merged.snps = 1188973;
config.frames.EUR_100K_1188K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 100000, fullfile(norstore_path, 'SYNGP/EUR_100K_1M_merged/all'));
config.frames.EUR_100K_1188K_merged.bim = PlinkRead_bim(fullfile(norstore_path, 'SYNGP/EUR_100K_1M_merged/all'));

% Matlab 1M reference template (w_hm3) -- 100K subjects
config.frames.EUR_100K_1190K_ref.subj = 100000;
config.frames.EUR_100K_1190K_ref.snps = 1190321;
config.frames.EUR_100K_1190K_ref.bim  = PlinkRead_bim(fullfile(norstore_path, 'SYNGP/1m_ref'), true, '%s %s %f %d %s %s %s %s');
