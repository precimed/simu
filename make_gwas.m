function frame = make_gwas(frame, config, varargin)
    % MAKE_GWAS generates gwas estimates of effect size and corresponding p-value

    p = inputParser;
    addOptional(p, 'snpstep', 100);    % how many SNPs to scan per one iteration
    addOptional(p, 'subj_idx', {});    % which subjects to involve in GWAS
    % subj_idx can be
    % - empty
    % - an array with subjects to include in gwas (will be used for all traits)
    % - an cell array of length ntraits; each element is the list to use for each trait
    
    parse(p, varargin{:}); opts = p.Results;
    
    if frame.subj ~= config.frames.(frame.name).subj, error('Number of subjects mismatch'); end;
    if frame.snps ~= config.frames.(frame.name).snps, error('Number of SNPs mismatch'); end;

    nsnp = frame.snps;
    nsubj = frame.subj;
    snpstep = opts.snpstep;
    ntraits = size(frame.phenotype, 2);
    reader = config.frames.(frame.name).reader;

    opts = fix_and_validate(frame, config, opts);
    
    gwasbeta = zeros(nsnp, ntraits);
    gwaspval = zeros(nsnp, ntraits);
    beginningOfTime = now;
    for i=1:snpstep:nsnp
        e = min(i+snpstep, nsnp);
        X = reader(1:nsubj, i:e); X = double(X); X = (X-repmat(mean(X,1),nsubj,1));
        % X = X ./ repmat(std(X), [nsubj, 1]);   ???????? do we need somehting like this or not?
        for trait=1:ntraits
            observedPhenotype = frame.phenotype(opts.subj_idx{trait}, trait);
            [rhovec_tmp, pvec_tmp] = corr(X(opts.subj_idx{trait}, :), double(observedPhenotype));
            if any(isnan(rhovec_tmp)) || any(isnan(pvec_tmp)), error('corr return nan\n'); end;
            gwasbeta(i:e, trait) = rhovec_tmp;
            gwaspval(i:e, trait) = pvec_tmp;
        end
        fprintf(1,'make_gwas: %.1f%% done, Now:%s eta:%s\n', 100 * (i+snpstep) / nsnp, datestr(now,'dddd HH:MM:SS'),datestr(beginningOfTime+(now-beginningOfTime)/(i+snpstep)*nsnp));
    end

    frame.gwasbeta = gwasbeta;
    frame.gwaspval = gwaspval;
    frame.gwassize = zeros(size(gwaspval)); for trait=1:ntraits, frame.gwassize(:, trait) = ones(frame.snps, 1) * length(opts.subj_idx{trait}); end
    frame.opts.make_gwas = opts;
end



function opts = fix_and_validate(frame, config, opts)
    nsubj = frame.subj;
    ntraits = size(frame.phenotype, 2);

    if isempty(opts.subj_idx),
        for trait=1:ntraits, opts.subj_idx{trait} = 1:nsubj; end;
    end;
    if ~iscell(opts.subj_idx),
        a=opts.subj_idx; opts.subj_idx=cell(1, ntraits);
        for trait=1:ntraits, opts.subj_idx{trait} = a; end;
    end
    if length(opts.subj_idx) ~= ntraits, error('Invalid argument: subj_idx'); end;
    for trait=1:ntraits
        if islogical(opts.subj_idx{trait}), opts.subj_idx{trait} = find(opts.subj_idx{trait}); end;
        if any(opts.subj_idx{trait} < 1 | opts.subj_idx{trait} > nsubj), error('Argument out of range: subj_idx'); end;
    end
end