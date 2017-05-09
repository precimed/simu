function frame = reframe(frame, config, framename)
    % REFRAME changes frame and align matrices to the new set of SNPs.

    oldcfg = config.frames.(frame.name);
    newcfg = config.frames.(framename);
    fn = fieldnames(frame);
    
    % Use information from bim files to rearrange SNPs
    if isfield(oldcfg, 'bim') && isfield(newcfg, 'bim') && (oldcfg.subj == newcfg.subj)
        % Changing the set of subjects is not allowed
        if length(oldcfg.bim.snplist) ~= frame.snps, error('Number of SNPs mismatch (in the old frame)'); end;
        if length(newcfg.bim.snplist) ~= newcfg.snps, error('Number of SNPs mismatch (in the new frame)'); end;
        [is_in_new_frame, index_to_new_frame] = ismember(oldcfg.bim.snplist, newcfg.bim.snplist);
        if ~all(is_in_new_frame), warning('Unable to map %i SNPs to the new frame; %i SNPs mapped successfully', sum(~is_in_new_frame), sum(is_in_new_frame)); end;
        for i=1:length(fn)
            if size(frame.(fn{i}), 1) ~= frame.snps, continue; end;
            frame.(fn{i}) = align_to_ref(frame.(fn{i}), length(newcfg.bim.snplist), index_to_new_frame);
        end
    end

    % Taking subset of subjects
    if (oldcfg.subj >= newcfg.subj)
        for i=1:length(fn)
            if size(frame.(fn{i}), 1) ~= frame.subj, continue; end;
            matrix = frame.(fn{i});
            matrix = matrix(1:newcfg.subj, :);
            frame.(fn{i}) = matrix;
        end
    end

    % Removing fields that will not match new frame
    to_be_removed = {};
    for i=1:length(fn)
        if size(frame.(fn{i}), 1) == frame.subj && oldcfg.subj ~= newcfg.subj, frame = rmfield(frame, fn{i}); to_be_removed{end+1} = fn{i}; end;
        if size(frame.(fn{i}), 1) == frame.snps && oldcfg.snps ~= newcfg.snps, frame = rmfield(frame, fn{i}); to_be_removed{end+1} = fn{i}; end;
    end
    if ~isempty(to_be_removed), warning('reframe is removing fields: %s\n', sprintf('%s ', to_be_removed{:})); end;

    frame.snps = newcfg.snps;
    frame.subj = newcfg.subj;
    frame.name = framename;
end

function aligned_vec = align_to_ref(vec, nsnp_ref, index_to_ref)
    numtraits = size(vec, 2);
    aligned_vec = nan(nsnp_ref, numtraits);
    aligned_vec(index_to_ref(index_to_ref ~= 0), :) = vec(index_to_ref ~= 0, :);
end
