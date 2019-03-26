% do power simulation for networks:
function [finalAUCs, samplesize] = PowerSim(m, s, maxn, startsize)
%
% do a plot of AUC values in hold-out cross validation for various values
% of N, the samplesize, for m predictors that have normal distributions
% dependent upon phenotype and that the difference in these is itself
% normally distributed with mean 0 and variance = s^2.
%
% BETA.
%
% Copyright 2013 Michael McGeachie.

if (nargin < 4)
    startsize = 0;
end

FOLDS = 3;
STEPSIZE = 50;
numreps = 5;
priorPrecision.nu = 25;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 25;
priorPrecision.maxParents = 2;

topfilter = 200;
dofiltering = true;

cols = cell(1,1+m);
cols{1} = 'pheno';
for i = 1:m
    cols{i+1} = ['normmet',num2str(i)];
end
disc = false(1,length(cols));
disc(1) = true;

for i = 1:(maxn-startsize)/STEPSIZE
    % generate dataset
    n = startsize+STEPSIZE*i;
    pheno = ones(n,1);
    pheno(ceil(n/2):end) = 0;

    fprintf(1,'Simulating N = %d\n',n);
    aucs = zeros(1,numreps);
    for j = 1:numreps
        % normally distributed data
        % should update this data to have variable standard deviations
        % ones that are similar to other observed metabolomics datasets
        data = randn(n,m);
        % differences conditional on phenotype are ~ N(0,s^2)
        diffs = randn(1,m) * s;
        data(pheno == 1,:) = repmat(diffs,sum(pheno ==1),1) + data(pheno == 1,:);
    
        cvdata = [pheno,data];
        if (dofiltering)
            % do a CV-network learning experiment
            [aucs(j), ~, ~] = CVModelLearnFilter(cvdata, cols, 'pheno', ...
                priorPrecision, FOLDS, 'simnormdata', true, true, 1, topfilter);
        else        
            % do a CV-network learning experiment
            [aucs(j), ~, ~] = CVModelLearnEx(cvdata, cols, 'pheno', ...
                priorPrecision, FOLDS, 'simnormdata', false, true, 1);
        end
    end
    
    % average the AUC for each of numreps replications
    finalAUCs{i} = aucs;
    samplesize(i) = n;
end

% plot/return the value of average AUC vs. n


