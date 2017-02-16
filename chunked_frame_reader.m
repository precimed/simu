function genomat = chunked_frame_reader(subjvec, snpvec, nsubj, snpdetailpath, genotypepath)
    % Read genotypes from a single file
    % subjvec --- indices of subjects to read
    % snpvec  --- indices of SNPs to read
    % nsubj   --- total number of subjects present in the file
    % bfile   --- path to plink BED file (without extension)
    %
    persistent PRS_chunk_file         % cell array [nchunks x 1], file name of each non-empty chunk
    persistent PRS_chunk_of_snp       % vector, [nsnp x 1],       which chunk SNP belongs to
    persistent PRS_chunk_index_of_snp % vector, [nsnp x 1],       unity-based index of each SNP in its chunk
    persistent PRS_snpdetailpath PRS_genotypepath
    
    % force rebuild persistent variables if snpdetailpath or genotypepath
    % has changed compared to the previous call of this function
    if isempty(PRS_snpdetailpath) || ~strcmp(PRS_snpdetailpath, snpdetailpath), PRS_chunk_file = []; end;
    if isempty(PRS_genotypepath)  || ~strcmp(PRS_genotypepath,  genotypepath),  PRS_chunk_file = []; end;
    
    if isempty(PRS_chunk_file)
        PRS_snpdetailpath = snpdetailpath; PRS_genotypepath = genotypepath;
        snp_detail = load(snpdetailpath);
        PRS_chunk_file = {};
        PRS_chunk_of_snp = nan(sum(snp_detail.snp_size(:)), 1);       
        PRS_chunk_index_of_snp = nan(sum(snp_detail.snp_size(:)), 1);  
        snp_index = 0;
        fprintf('Verifying integrity of chunks... ');
        for i=1:snp_detail.nchrom
        for j=1:snp_detail.nchunk(i)
            fileprefix=sprintf('%s/chr%d_chunk%d',genotypepath,i,j);
            bedprefix=[fileprefix,'.bed'];
            if ~exist(bedprefix,'file') && (snp_detail.snp_size(i, j) > 0)
                error('123');
            end
            snps_per_chunk = snp_detail.snp_size(i, j);
            if snp_detail.snp_size(i, j) > 0
                if ~exist(bedprefix,'file'), error('Expected %s file not found --- %i SNPs lost', bedprefix, snps_per_chunk); end
                PRS_chunk_file{end+1, 1} = fileprefix;
                PRS_chunk_of_snp((snp_index + 1) : (snp_index + snps_per_chunk)) = length(PRS_chunk_file);
                PRS_chunk_index_of_snp((snp_index + 1) : (snp_index + snps_per_chunk)) = 1:snps_per_chunk;
                snp_index = snp_index + snps_per_chunk;
            end
        end
        end
        fprintf('OK. Found %i chunks with %i SNPs in total.', length(PRS_chunk_file), snp_index);
    end
    
    if islogical(snpvec), snpvec = find(snpvec); end;
    if ~issorted(snpvec), error('chunked_frame_reader expect sorted list of snps'); end;

    genomat = zeros(length(subjvec),length(snpvec),'int8');
    chunks_to_scan = unique(PRS_chunk_of_snp(snpvec));
    for chunk = chunks_to_scan'
        snpvec_for_this_chunk = snpvec(PRS_chunk_of_snp(snpvec) == chunk);
        tmpgeno = PlinkRead_binary2(nsubj, PRS_chunk_index_of_snp(snpvec_for_this_chunk), PRS_chunk_file{chunk});
        tmpgeno = tmpgeno(subjvec, :);
        genomat(:, PRS_chunk_of_snp(snpvec) == chunk) = tmpgeno;
    end
end
