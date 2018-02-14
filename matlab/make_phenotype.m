function frame = make_phenotype(frame, config, varargin)
    % MAKE_PHENOTYPE adjust truepheno to a given level of heritability

    p = inputParser;
    addOptional(p, 'h2', 0.8);         % how many SNPs to scan per one iteration
    parse(p, varargin{:}); opts = p.Results;

    nsubj = frame.subj;
    h2 = opts.h2;
    ntraits = size(frame.truepheno, 2);
    
    for trait = 1:ntraits
        truePhenotype = frame.truepheno(:, trait);
        noisevec = randn(nsubj, 1);
        if h2 == 0 || var(truePhenotype) == 0
          observedPhenotype = noisevec / sqrt(var(noisevec));
          % This is an interesting case --- why there are some true casual SNPs, and yet heritability is 0?
          % This imply that environment has far larger effect that those SNPs.
          % Mathematically in this case the factor (1-h2)/h2 would go to infinity,
          % So we just pretend that observed phenotype has no relation with true phenotype (thus with beta vec).
        else
          observedPhenotype = truePhenotype + noisevec*sqrt((1-h2)/h2*var(truePhenotype)/var(noisevec));
          %=> var(truePhenotype) / var(observedPhenotype) == h2
        end
        frame.phenotype(:, trait) = observedPhenotype;
    end
    
    frame.opts.make_phenotype = opts;
end