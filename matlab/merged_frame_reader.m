function genomat = merged_frame_reader(subjvec, snpvec, nsubj, bfile)
    % Read genotypes from a single file
    % subjvec --- indices of subjects to read
    % snpvec  --- indices of SNPs to read
    % nsubj   --- total number of subjects present in the file
    % bfile   --- path to plink BED file (without extension)
    
    genomat = PlinkRead_binary2(nsubj, snpvec, bfile);  % int8, nsubj * nsnp
    genomat = genomat(subjvec, :);
end
