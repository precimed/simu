function frame = make_truebeta_gmm(frame, config, varargin)
    % MAKE_TRUEBETA_GMM generates two vectors of true effect sizes from gaussian mixture model
    % Example:
    % frame.snps = 1000; make_truebeta_gmm(frame, 'pivec', [0.01], 'sig1vec', [1.0], 'sig2vec', [1.0], 'rhovec', [0.5])

    p = inputParser;
    addOptional(p, 'pivec', []);    % proportion of snps in each component (except null)
    addOptional(p, 'sig1vec', []);  % variance (first component)
    addOptional(p, 'sig2vec', []);  % variance (second component)
    addOptional(p, 'rhovec', []);   % correlation (between two components)
    addOptional(p, 'mask', []);     % restrict heritability to certain SNPs
                                    % mask can be either logical array of size "snps x 1"
                                    % or an array of indices between 1 and snps.
    parse(p, varargin{:}); opts = p.Results;

    mask = opts.mask;
    pivec = opts.pivec;
    sig1vec = opts.sig1vec;
    sig2vec = opts.sig2vec;
    rhovec = opts.rhovec;

    nsnp = frame.snps;
    if isempty(mask), mask = 1:nsnp; end;
    if islogical(mask), mask = find(mask); end;

    if nsnp <= 0, error('nsnp must be positive'); end;
    if sum(floor(pivec * nsnp)) > nsnp, error('sum of pivec must be less than 1'); end;
    if sum(floor(pivec * nsnp)) > length(mask), error('sum of pivec is to large for given mask'); end;
    if length(unique([length(sig1vec), length(sig2vec), length(rhovec), length(pivec)])) > 1,
        error('pivec, sig1vec, sig2vec, rhovec must have the same length');
    end;

    betamat = zeros(nsnp, 2);
    mix = zeros(nsnp, 1);
    num_mixtures = length(pivec);

    idx = mask(randperm(length(mask))); idx_start = 1;
    for i=1:num_mixtures
        mu = [0 0];
        sigma = [sig1vec(i), sqrt(sig1vec(i) * sig2vec(i)) * rhovec(i);
                 sqrt(sig1vec(i) * sig2vec(i)) * rhovec(i), sig2vec(i)];
        num_causals = floor(nsnp * pivec(i));
        %idx = randsample(1:nsnp, num_causals); betamat(idx, :) = betamat(idx, :) + mvnrnd(mu, sigma, num_causals);

        idx_end = idx_start + num_causals - 1;
        betamat(idx(idx_start:idx_end), :) = mvnrnd(mu, sigma, num_causals);
        mix(idx(idx_start:idx_end), :) = i;
        idx_start = idx_start + num_causals;
    end
    
    frame.opts.make_truebeta = opts;
    frame.truebeta      = betamat;
    frame.mixture       = repmat(mix, [1 2]);
end





