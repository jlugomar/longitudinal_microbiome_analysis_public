function BFs = BFPriorExplore(data, cols, pheno, nvals, disc)
%BFs = BFPriorExplore(data, cols, pheno, nvals, disc)
% 
% Test impact of priors on bayes factors for a given dataset
%
% This computes bayes factors for each column of non-phenotype data in the
% input, using an increasingly strong prior.  It starts with 
%   priorPrecision.alpha = 1
%   priorPrecision.nu = 1
% And computes BFs using these numbers.  For NVALS interations it increases
% these priorPreicsion values and computes further BFs.  Plots log BFs for each
% column on a graph.
%
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to compute Bayes Factors 
%   against.  Is matched against the COLS array.
% NVALS: number of different values of priorPrecision.nu and
%   priorPreicions.alpha to simulate.  Going from [1 to NVALS].
%
% OUTPUT:
% BFs: array of log bayes factors matching the cols.
%
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.


if (nargin < 4)
    nvals = 50;
end
if (nargin < 5)
    disc = IsDiscrete(data);
end

priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 1;
priorPrecision.maxParents = 2;

BFs = zeros(nvals,length(cols));
for i = 1:nvals
    priorPrecision.nu = i * 1;
    priorPrecision.alpha = i * 1;
    ret = BayesFactorScore(data, cols, pheno, priorPrecision, disc);
    BFs(i,:) = ret';
end

% make a figure
figure();
maxcols = size(data,2);
cc = colormap(jet(maxcols));
hold on;
for i = 1:maxcols
    plot(BFs(:,i),'color',cc(i,:));
    [m,mind] = max(BFs(:,i));
    plot(mind,m,'o','color',cc(i,:));
end
hold off;
title('BFs vs. Prior Strength');
ylabel('Log Bayes Factor');
xlabel('prior virtual sample size');















