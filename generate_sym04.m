% input:
% pi1u, pi2u - total polygenicity of each trait
% pi12 - polygenic overlap (pi12 <= pi1, pi12 <= pi2)
% h2 - heritability
% gencorr - genetic correlation in the pleiotropic component
% subjoverlap - 'zero', 'partial', 'full'
%
% example:
% pi1u=1e-5; pi2u=1e-4; pi12=1e-6; h2=0.5; gencorr=0.6; subjoverlap='zero';
% filename_id={'QRDKXL'; 'ASCPL'};

pi1 = pi1u - pi12;
pi2 = pi2u - pi12;
if (pi1 < 0 || pi2 < 0 || pi12 < 0 || (pi1+pi2+pi12 > 1))
    error('error in pi1u, pi2u, or pi12');
end

if (h2 < 0 || h2 > 1), error('error in h2'); end;
if (gencorr < -1 || gencorr > 1), error('error in gencorr'); end;

if ~exist('config', 'var')
    t = 'H:\NORSTORE\SYNGP';                       if exist(t, 'dir'), norstore_path = t; end;
    t = '/usit/abel/u1/oleksanf';                  if exist(t, 'dir'), norstore_path = t; end;
    t = 'E:\EUR_100K_9M_merged';                   if exist(t, 'dir'), EUR_100K_9M_merged_path = t; end;
    t = '/work/users/oleksanf/EUR_100K_9M_merged'; if exist(t, 'dir'), EUR_100K_9M_merged_path = t; end;

    % A subset of 80M and 9M template
    config.frames.EUR_100K_8801K_merged.subj = 100000;
    config.frames.EUR_100K_8801K_merged.snps = 8801249;
    config.frames.EUR_100K_8801K_merged.reader = @(subjvec, snpvec)merged_frame_reader(subjvec, snpvec, 100000, fullfile(EUR_100K_9M_merged_path, 'all'));
    config.frames.EUR_100K_8801K_merged.bim = PlinkRead_bim(fullfile(EUR_100K_9M_merged_path, 'all'));
    config.frames.EUR_100K_8801K_merged.bim = rmfield(config.frames.EUR_100K_8801K_merged.bim, {'chrvec', 'cMvec', 'bpvec', 'A1vec', 'A2vec'});

    % Matlab 9M template
    config.frames.EUR_100K_9279K_ref.subj = 100000;
    config.frames.EUR_100K_9279K_ref.snps = 9279485;
    config.frames.EUR_100K_9279K_ref.bim  = PlinkRead_bim(fullfile(norstore_path, '9279485_ref'), true, '%s %s %f %d %s %s %s %s');
    config.frames.EUR_100K_9279K_ref.bim = rmfield(config.frames.EUR_100K_9279K_ref.bim, {'chrvec', 'cMvec', 'bpvec', 'A1vec', 'A2vec'});

    % Matlab 1M reference template (w_hm3)
    config.frames.EUR_100K_1190K_ref.subj = 100000;
    config.frames.EUR_100K_1190K_ref.snps = 1190321;
    config.frames.EUR_100K_1190K_ref.bim  = PlinkRead_bim(fullfile(norstore_path, '1m_ref'), true, '%s %s %f %d %s %s %s %s');
end
fieldnames(config.frames)

for repeat_index = 1:length(filename_id)
rng('shuffle')

frame = make_empty_frame (config, 'EUR_100K_8801K_merged');

nsubj = frame.subj;
s500 = floor(nsubj/2); s501 = s500+1; s1000 = nsubj;
s333 = floor(nsubj/3); s334 = s333+1; s667 = 2*s333+1;
overlap=[];
overlap.zero = {1:s500, s501:s1000}; % zero overlap
overlap.partial = {1:s667, s334:s1000}; % partial overlap
overlap.full = {1:s1000, 1:s1000};   % full overlap

overlap = overlap.(subjoverlap);

pivec = []; sig1vec = []; sig2vec = []; rhovec = [];
if (pi1 > 0) pivec = [pivec pi1]; sig1vec = [sig1vec 1]; sig2vec = [sig2vec 0]; rhovec = [rhovec 0]; end;
if (pi2 > 0) pivec = [pivec pi2]; sig1vec = [sig1vec 0]; sig2vec = [sig2vec 1]; rhovec = [rhovec 0]; end;
if (pi12 > 0) pivec = [pivec pi12]; sig1vec = [sig1vec 1]; sig2vec = [sig2vec 1]; rhovec = [rhovec gencorr]; end;

frame = make_truebeta_gmm(frame, config, 'pivec', pivec, 'sig1vec', sig1vec, 'sig2vec', sig2vec, 'rhovec', rhovec);
frame = make_truepheno   (frame, config, 'snpstep', 100);

frame = make_phenotype(frame, config, 'h2', h2);
frame = make_gwas     (frame, config, 'subj_idx', overlap);
frame = make_gwaslogp (frame, config);

templates = {'EUR_100K_9279K_ref', 'EUR_100K_1190K_ref'};
template_short_name  = {'9M', '1M'};
for it = 1:length(templates)
    frame1 = reframe(frame, config, templates{it});
    filename = sprintf('template=%s_pi1u=%.0e_pi2u=%.0e_pi12=%.0e_gencorr=%.2f_h2=%.2f_subjoverlap=%s_id=%s.mat', template_short_name{it}, pi1u, pi2u, pi12, gencorr, h2, subjoverlap, filename_id{repeat_index});
    save(filename, '-struct', 'frame1', 'truebeta', 'mixture', 'truepheno', 'phenotype', 'gwasbeta', 'gwaspval', 'nvec', 'logpvec', 'zvec', '-v7.3');
end

end