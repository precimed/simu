function plot_qq(logpvec, varargin)
    % PLOT_QQ generates standard QQ plot

    p = inputParser;
    addOptional(p, 'qqlim', 7);
    addOptional(p, 'maxpoints', 500000);
    parse(p, varargin{:}); opts = p.Results;

    logpvec = logpvec(~isnan(logpvec));
    n=length(logpvec); 
    y = sort(logpvec);
    x=-log10(1-(1:n)/n)';
    
    mp = opts.maxpoints;  % limit number of display points to avoid matlab plotting issues
    if length(x) > mp, x=x(end-mp:end); y=y(end-mp:end); end;
    
    hold on
    plot(x, y, 'b');
    plot([0 opts.qqlim], [0 opts.qqlim], 'k--');
    axis([0 opts.qqlim 0 opts.qqlim]);
end
