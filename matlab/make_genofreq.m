function frame = make_genofreq(frame, config, varargin)
    % MAKE_GENOFREQ calculate frequencies of AA, Aa, and aa genotypes.

    p = inputParser;
    addOptional(p, 'snpstep', 100);  % how many SNPs to scan per one iteration
    parse(p, varargin{:}); opts = p.Results;

    if frame.subj ~= config.frames.(frame.name).subj, error('Number of subjects mismatch'); end;
    if frame.snps ~= config.frames.(frame.name).snps, error('Number of SNPs mismatch'); end;

    nsnp = frame.snps;
    nsubj = frame.subj;
    reader = config.frames.(frame.name).reader;
    snpstep = opts.snpstep;

    genofreq = zeros(nsnp, 3);
    beginningOfTime=now;
    for i=1:snpstep:nsnp
        e = min(i+snpstep, nsnp);
        X = reader(1:nsubj, i:e);
        genofreq(i:e, 1) = sum(X == 0);
        genofreq(i:e, 2) = sum(X == 1);
        genofreq(i:e, 3) = sum(X == 2);
        fprintf(1,'make_genofreq: %.1f%% done, Now:%s eta:%s\n', 100 * (i+snpstep) / nsnp, datestr(now,'dddd HH:MM:SS'),datestr(beginningOfTime+(now-beginningOfTime)/(i+snpstep)*nsnp));
    end

    frame.genofreq = genofreq;
    frame.opts.make_genofreq = opts;
end
