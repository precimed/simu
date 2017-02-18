function frame = make_truepheno(frame, config, varargin)
    % MAKE_TRUEPHENO generates true phenotypes (without adding environmental noise)

    p = inputParser;
    addOptional(p, 'snpstep', 100);  % how many SNPs to scan per one iteration
    addOptional(p, 'stdnorm', false);  % normalize genotypes as in LD score regression (=> rare SNPs has larger effect size)
    parse(p, varargin{:}); opts = p.Results;

    if frame.subj ~= config.frames.(frame.name).subj, error('Number of subjects mismatch'); end;
    if frame.snps ~= config.frames.(frame.name).snps, error('Number of SNPs mismatch'); end;

    nsubj = frame.subj;
    snpstep = opts.snpstep;
    ntraits = size(frame.truebeta, 2);
    reader = config.frames.(frame.name).reader;
    
    frame.truepheno = nan(nsubj, ntraits);
    beginningOfTime = now;
    for trait = 1:ntraits
        betavec = frame.truebeta(:, trait);
        betavec_idx = find(betavec ~= 0);
        truePhenotype = zeros(nsubj, 1);

        for i=1:snpstep:length(betavec_idx)
            e = min(i+snpstep, length(betavec_idx));
            genotypes = reader(1:nsubj, betavec_idx(i:e));
            X = double(genotypes); X = (X-repmat(mean(X,1),nsubj,1));
            if opts.stdnorm, X = X ./ repmat(std(X), [nsubj, 1]); end;
            truePhenotype = truePhenotype + X * betavec(betavec_idx(i:e));
            fprintf(1,'make_truepheno (trait %i of %i): %.1f%% done, Now:%s eta:%s\n', trait, ntraits, 100 * (i+snpstep) / length(betavec_idx), datestr(now,'dddd HH:MM:SS'),datestr(beginningOfTime+(now-beginningOfTime)/(i+snpstep)*length(betavec_idx)));
        end

        frame.truepheno(:, trait) = truePhenotype;
    end
    
    frame.opts.make_truepheno = opts;
end