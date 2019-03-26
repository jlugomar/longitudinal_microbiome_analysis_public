function y = PCASNPs(snps, kval)
% y = PCASNPs(snps, kval)
%
% Computes principal components of a matrix of SNPs.  Returns the loadings
% of the cohort onto the top PCs, where there are enough PCs to account for
% a KVAL-fraction of the variance in the SNP matrix.
%
% INPUT:
% SNPS: matrix of SNPs
% KVAL: value in [0,1], representing how many of the top PCs to return.
%
% OUTPUT:
% Y: loadings of each row of SNP data onto the top PCs.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

% center the SNP data:
csnps = single(snps) - repmat(mean(snps),size(snps,1),1);

% first calculate covariance matrix
c = cov(csnps);

% find the top NP pca's of c:
[v,d] = eig(c);

% sort the PCA's 
[s,sorder] = sort(diag(d),'descend');
%d = d(sorder, sorder);
v = v(:, sorder);

% compute the "energy distribution" of the eigen values
g = cumsum(s);
n = find(g > kval * max(g),1);
w = v(:,1:n);

% standard deviations
sds = sqrt(diag(c));
divisor = repmat(sds',size(csnps,1),1);
zsnps = zeros(size(csnps));
zsnps(divisor ~= 0) = csnps(divisor ~= 0) ./ divisor(divisor ~= 0);

y = w' * zsnps';