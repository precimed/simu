function betamat = make_gmm_beta(nsnp, pivec, sig1vec, sig2vec, rhovec)
    % MAKE_GMM_PHENO generates two vectors of effect sizes from gaussian mixture model
    % make_gmm_beta(10000, [0.01], [1.0], [1.0], [0.5])

    if nsnp <= 0, error('nsnp must be positive'); end;
    if sum(floor(pivec * nsnp) > nsnp), error('sum of pivec must be less than 1'); end;
    if length(unique([length(sig1vec), length(sig2vec), length(rhovec), length(pivec)])) > 1,
        error('pivec, sig1vec, sig2vec, rhovec must have the same length');
    end;

    betamat = zeros(nsnp, 2);
    num_mixtures = length(pivec);

    idx = randperm(nsnp); idx_start = 1;
    for i=1:num_mixtures
        mu = [0 0];
        sigma = [sig1vec(i), sqrt(sig1vec(i) * sig2vec(i)) * rhovec(i);
                 sqrt(sig1vec(i) * sig2vec(i)) * rhovec(i), sig2vec(i)];
        num_causals = floor(nsnp * pivec(i));
        %idx = randsample(1:nsnp, num_causals); betamat(idx, :) = betamat(idx, :) + mvnrnd(mu, sigma, num_causals);

        idx_end = idx_start + num_causals - 1;
        betamat(idx(idx_start:idx_end), :) = mvnrnd(mu, sigma, num_causals);
        idx_start = idx_start + num_causals;
    end
end
