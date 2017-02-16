function frame = make_empty_frame(config, framename)
    frame = struct();
    frame.snps = config.frames.(framename).snps;
    frame.subj = config.frames.(framename).subj;
    frame.name = framename;
end