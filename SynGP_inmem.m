function [betavec, rhovec, pvec] = SynGP_inmem(G, nsubj, nsnp, pi1, h2, Nfrac, snpstep, betavec)

if ~exist('snpstep', 'var'), snpstep = 10000; end;

if ~exist('betavec', 'var'),
    nsnpCausal = max(floor(pi1*nsnp), 1);
    if nsnpCausal > 1
        if ~mod(nsnpCausal,2), nsnpCausal = nsnpCausal + 1;  end
        nsnpCausal = min(nsnpCausal, nsnp);
        cdfUniformVec = linspace(0,1,nsnpCausal+2);
        cdfUniformVec = cdfUniformVec(2:end-1);
        betaVecCausals = norminv(cdfUniformVec);
        betaVecCausalsNorm = betaVecCausals/std(betaVecCausals);  % ==> var(betaVecCausalsNorm)==1
    else
        betaVecCausalsNorm = [ 1 ];
    end
    permvec = randperm(nsnp);
    betavec = zeros(nsnp, 1);
    betavec(permvec(1:nsnpCausal),1) = betaVecCausalsNorm;
end
betavec_idx = find(betavec ~= 0);

truePhenotype = zeros(nsubj, 1);

for i=1:snpstep:length(betavec_idx)
    fprintf('.');
    e = min(i+snpstep, length(betavec_idx));
    X = double(G(betavec_idx(i:e))); X = (X-repmat(mean(X,1),nsubj,1));
    truePhenotype = truePhenotype + X * betavec(betavec_idx(i:e));
end

permvec = randperm(nsubj);
indvecs = permvec(1:round(Nfrac*nsubj));

noisevec = randn(nsubj, 1);
if h2 == 0
  observedPhenotype = noisevec / sqrt(var(noisevec));
  % This is an interesting case --- why there are some true casual SNPs, and yet heritability is 0?
  % This imply that environment has far larger effect that those SNPs.
  % Mathematically in this case the factor (1-h2)/h2 would go to infinity,
  % So we just pretend that observed phenotype has no relation with true phenotype (thus with beta vec).
else
  observedPhenotype = truePhenotype + noisevec*sqrt((1-h2)/h2*var(truePhenotype)/var(noisevec));
  %=> var(truePhenotype) / var(observedPhenotype) == h2
end

rhovec = zeros(nsnp, 1);
pvec   = zeros(nsnp, 1);
for i=1:snpstep:nsnp
    fprintf('.');
    e = min(i+snpstep, nsnp);
    X = double(G(i:e)); X = (X-repmat(mean(X,1),nsubj,1));
    [rhovec_tmp, pvec_tmp] = corr(X(indvecs, :), double(observedPhenotype(indvecs)));
    rhovec(i:e) = rhovec_tmp;
    pvec(i:e) = pvec_tmp;
end

fprintf('OK. \n');

end

