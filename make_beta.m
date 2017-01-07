function betavec = make_beta(nsnp, pi1)
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
