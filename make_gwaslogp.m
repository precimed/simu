function frame = make_gwaslogp(frame, config, varargin)
    % MAKE_LOGPVEC generates standard columns such as logpvec and zvec

    p = inputParser;
    parse(p, varargin{:}); opts = p.Results;
    
    gwaslogp = -log10(frame.gwaspval);
    gwaslogp(isinf(gwaslogp)) = -log10(realmin);
    gwaszscore = -norminv(10.^-gwaslogp/2).*sign(frame.gwasbeta);

    frame.gwaslogp = gwaslogp;
    frame.gwaszscore = gwaszscore;
    frame.opts.make_gwaslogp = opts;
end