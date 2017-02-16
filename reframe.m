function frame = reframe(frame, config, framename)
    % REFRAME changes frame and align matrices to the new set of SNPs.

    oldcfg = config.frames.(frame.name);
    newcfg = config.frames.(framename);
    
    % Use information from bim files
    if isfield(oldcfg, 'bim') && isfield(newcfg, 'bim')
        if length(oldcfg.bim.snplist) ~= frame.snps, error('Number of SNPs mismatch'); end;
        [is_in_new_frame, index_to_new_frame] = ismember(oldcfg.bim.snplist, newcfg.bim.snplist);
        if ~all(is_in_new_frame), error('Unable to map some SNPs to the new frame'); end;
        fn = fieldnames(frame);
        for i=1:length(fn)
            if size(frame.(fn{i}), 1) ~= frame.snps, continue; end;
            frame.(fn{i}) = align_to_ref(frame.(fn{i}), length(newcfg.bim.snplist), index_to_new_frame);
        end
        frame.snps = length(newcfg.bim.snplist);
        frame.subj = newcfg.subj;
        frame.name = framename;
    else
        error('Unable to reframe from %s to %s', frame.name, framename);
    end
end

function aligned_vec = align_to_ref(vec, nsnp_ref, index_to_ref)
    numtraits = size(vec, 2);
    aligned_vec = nan(nsnp_ref, numtraits);
    aligned_vec(index_to_ref, :) = vec;
end
