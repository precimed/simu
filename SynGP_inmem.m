function [zvec, nvec, betavec, rhovec, logpvec, mix] = SynGP_inmem(bfile, reffile, nsubj, nsnp, betavec, h2, subj_idx, snpstep, mix, matfile)
% SynGP_inmem converts true effect size (betavec) into GWAS estimate (rhovec, pvec)
% bfile --- plink b-file to read genotypes, example: 'EUR_10K_2M_merged/all'
% reffile --- superset of SNPs to align the result

persistent prev_bfile prev_reffile index_to_ref nsnp_ref
bfile_changed = isempty(prev_bfile) || ~strcmp(prev_bfile, bfile);
reffile_changed = isempty(prev_reffile) || ~strcmp(prev_reffile, reffile);

if bfile_changed || reffile_changed
    % 2.5M template: align to 2558411 reference template; 2518104 SNPs overlap (98.42%)
    % 1M template:   align to 1190321 reference template; 1188973 SNPs overlap (99.89%)
    bim = PlinkRead_bim(bfile);
    if nsnp ~= length(bim.snplist), error('expected nsnp : %i', length(bim.snplist)); end;

    bim_ref = PlinkRead_bim(reffile, true, '%s %s %f %d %s %s %s %s');
    nsnp_ref = length(bim_ref.snplist);

    [in_ref, index_to_ref] = ismember(bim.snplist, bim_ref.snplist);
    %[in_index, ref_to_index] = ismember(bim_ref.snplist, bim.snplist);
    if ~all(in_ref), error('Unknown SNPs detected in the data'); end;
    prev_bfile = bfile;
    prev_reffile = reffile;
    
    clear bim bim_ref in_ref
end

if ~exist('betavec', 'var'), error('betavec is required'); end;
if ~exist('snpstep', 'var'), snpstep = 10000; end;
if ~exist('subj_idx', 'var'), subj_idx = 1:nsubj; end;
if any(size(betavec) ~= [nsnp, 1]), error('betavec shape should be #snps, 1'); end;
if islogical(subj_idx), subj_idx = find(subj_idx); end;
if any(subj_idx <= 0 | subj_idx > nsubj), error('subj_idx out of range'); end;

G = @(snpvec)PlinkRead_binary2(nsubj,snpvec,bfile);  % int8, nsubj * nsnp

betavec_idx = find(betavec ~= 0);
truePhenotype = zeros(nsubj, 1);

for i=1:snpstep:length(betavec_idx)
    %fprintf('.');
    e = min(i+snpstep, length(betavec_idx));
    X = double(G(betavec_idx(i:e))); X = (X-repmat(mean(X,1),nsubj,1));
    X = X ./ repmat(std(X), [nsubj, 1]);
    truePhenotype = truePhenotype + X * betavec(betavec_idx(i:e));
end

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

rhovec = zeros(nsnp, 1);
pvec   = zeros(nsnp, 1);
for i=1:snpstep:nsnp
    %fprintf('.');
    e = min(i+snpstep, nsnp);
    X = double(G(i:e)); X = (X-repmat(mean(X,1),nsubj,1));
    X = X ./ repmat(std(X), [nsubj, 1]);
    [rhovec_tmp, pvec_tmp] = corr(X(subj_idx, :), double(observedPhenotype(subj_idx)));
    if any(isnan(rhovec_tmp)) || any(isnan(pvec_tmp)), error('corr return nan\n'); end;
    rhovec(i:e) = rhovec_tmp;
    pvec(i:e) = pvec_tmp;
end

fprintf('OK. \n');

betavec = align_to_ref(betavec, nsnp_ref, index_to_ref);
rhovec  = align_to_ref(rhovec, nsnp_ref, index_to_ref);
pvec    = align_to_ref(pvec, nsnp_ref, index_to_ref);
mix     = align_to_ref(mix, nsnp_ref, index_to_ref);

logpvec = -log10(pvec);
logpvec(isinf(logpvec)) = -log10(realmin);
zvec = -norminv(10.^-logpvec/2).*sign(rhovec);
nvec = ones(size(zvec)) * length(subj_idx);

if exist('matfile', 'var')
    save(matfile, 'betavec', 'rhovec', 'logpvec', 'zvec', 'nvec', 'h2', '-v7.3');
    fprintf('Results are saved to %s\n', matfile);
end

end

function aligned_vec = align_to_ref(vec, nsnp_ref, index_to_ref)
    if length(vec) ~= length(index_to_ref), error('internal error'); end;
    aligned_vec = nan(nsnp_ref, 1);
    aligned_vec(index_to_ref) = vec;
end

function bim = PlinkRead_bim(fileprefix, header, format)

if ~exist('format', 'var'), format = '%s %s %f %d %s %s'; end;
if ~exist('header', 'var'), header = false; end;

bimprefix = [fileprefix,'.bim'];
fprintf('Read in plink bim from %s \r\n', bimprefix);
bimid = fopen(bimprefix,'r');
if header, fgets(bimid); end;
bimdata = textscan(bimid, format);
fclose(bimid);
bim.chrvec = bimdata{1};
bim.snplist = bimdata{2};
bim.cMvec = bimdata{3};
bim.bpvec = bimdata{4};
bim.A1vec = bimdata{5};
bim.A2vec = bimdata{6};

end

function genomat = PlinkRead_binary2(nsubj,snps,fileprefix)
persistent geno_values

% Written by Chun 2015
if ~issorted(snps), error('PlinkRead_binary2 expect sorted list of snps'); end;
nsnp = length(snps);

% bit shift to generate the genovalue matrix
bedprefix = [fileprefix,'.bed'];

if isempty(geno_values)
    geno_values = zeros(256,4,'int8');
    geno_code = [-1,1,0,2];
    shiftind = [0,2,4,6];
    indvec=zeros(1,4);

    for i = 1:256;
        ishift = int16(i-1);
        for j = 1:4;
            indvec(j) = bitand(bitsra(ishift,shiftind(j)),3) ;
        end
        indvec(indvec == 0) = 4;
        geno_values(i,:) = geno_code(indvec);
    end
end

% Read in the binary file
bedid = fopen(bedprefix);
genobin = uint16(fread(bedid, 3));

% Check magic number
if genobin(1) ~= 108;
	error('- Not a valid Plink BED file \r\n');
elseif genobin(2) ~= 27;
	error('- Not a valid Plink BED file \r\n');
elseif genobin(3) ~= 1;
	error('- Not in SNP-major format \r\n');
end

n_bytes = ceil(nsubj/4);
genomat = zeros(nsubj,nsnp,'int8');
for i = 1:nsnp;
    fseek(bedid, 3 + (snps(i) - 1) * n_bytes, 'bof');
    genobin = uint16(fread(bedid, n_bytes));
    if length(genobin) ~= n_bytes, error('-- Invalid number of entries from the bed \r\n'); end
    tmp_values = geno_values(genobin + 1, :)';
    tmp_values = tmp_values(:);
    genomat(:,i) = tmp_values(1:nsubj);
end
fclose(bedid);

end
