function [bestParams,aucs] = CVExploreParams(data, cols, pheno, experimentname)
% repeatedly calls CVModelLearn() with differen parameter settings, then
% compares those settings on AUC; returning the parameters that obtain the
% greatest AUC in cross-validation.
%
% among the params to optimize:
%   priorPrecision.nu; % prior sample size for prior variance estimate
%   priorPrecision.sigma2; % prior variance estimate
%   priorPrecision.alpha; % prior sample size for discrete nodes
%   priorPrecision.maxParents; % hard-limit on the number of parents each
%       node can have
%
% copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.

FOLDS = 5;
VERBOSE = false;

% loop through all these values: 
%nuVals = [5, 10, 25, 50, 100];
%sigmaVals = [.3, .5, 1, 2, 3.5, 5, 10];
%alphaVals = [10, 20, 50, 100, 250];
%maxparentsVals = [2,3];

nuVals = [5, 10, 25, 50, 100];
sigmaVals = [.5, 1, 2, 3.5, 5, 10];
alphaVals = [10, 20, 50, 100, 250];
maxparentsVals = [2];

% testing vals
%nuVals = [1, 5];
%sigmaVals = [.1, .3];
%alphaVals = [10, 20];
%maxparentsVals = [2, 3];

n = length(nuVals) * length(sigmaVals) * length(alphaVals) * length(maxparentsVals);
aucs = zeros(1,n);
inds = zeros(1,4);
maxVals = [length(nuVals), length(sigmaVals), length(alphaVals), length(maxparentsVals)];
for i = 1:n
    % this uses inds + 1 because the bitincbasevec() function counts in
    % arbitrary base starting with 0 in each digit.  So just add + 1 here
    % to get the right array indices
    priorPrecision.nu = nuVals(inds(1)+1);
	priorPrecision.sigma2 = sigmaVals(inds(2)+1);
    priorPrecision.alpha = alphaVals(inds(3)+1);
    priorPrecision.maxParents = maxparentsVals(inds(4)+1);
    fprintf(1,'Model %d has: Nu:%d sigma:%f alpha:%d maxParents:%d \n',...
        i, priorPrecision.nu, priorPrecision.sigma2, priorPrecision.alpha, priorPrecision.maxParents);
    [aucs(i), ~, ~] = CVModelLearnEx(data, cols, pheno, priorPrecision, FOLDS, ...
        experimentname, VERBOSE, false, true, false);
    fprintf(1,'\t Has AUC: %f\n',aucs(i));
    inds = bitincbasevec(inds, maxVals);
    if (isempty(inds))
        break;
    end
end

% find maximum:
[~,mind] = max(aucs);
maxInds = zeros(1,4);
for i = 1:mind-1
    maxInds = bitincbasevec(maxInds, maxVals);
end
bestParams.nu = nuVals(maxInds(1)+1);
bestParams.sigma2 = sigmaVals(maxInds(2)+1);
bestParams.alpha = alphaVals(maxInds(3)+1);
bestParams.maxParents = maxparentsVals(maxInds(4)+1);



