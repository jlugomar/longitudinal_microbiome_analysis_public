function BFs = BayesFactorScore(data, cols, pheno, priorPrecision, disc, verbose)
%bfs = BayesFactorScore(data, cols, pheno, priorPrecision)
%
% This function computes Bayes Factors for the dependence of each variable
% with the phenotype variable.  This is suitable for filtering for domains 
% with too many variables to be considered by usual Bayes Network methods. 
%
% The Bayes Factor is a bayesian likelihood ratio test that computes the
% ratio of posterior probabilites of two quantities: 1) the probability of
% the variable being statistically dependent upon the phenotype, and 2) the
% probability of that variable being independent of the phenotype.  Given
% in log scale.  Values > 0 mean the variable is more likely to be
% associated than not.
% 
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to compute Bayes Factors 
%   against.  Is matched against the COLS array.
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%
% OUTPUT:
% BFS: array of log bayes factors matching the cols.
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    verbose = true;
end
if (nargin < 4)
    priorPrecision.alpha = 10; %prior frquency for discrete data
    priorPrecision.nu = 10; %prior sample size for continuous data
    priorPrecision.sigma2 = 1; %variance in the prior (hypothetical) samples for continuous data
end

% find phenotype column
phencol = find(strcmp(pheno, cols));
pheno = data(:,phencol);

if (nargin < 5)
    % check for discrete data
    MAXDISCVALS = 4;
    disc = IsDiscrete(data,MAXDISCVALS);
end
n = size(data,2);

% compute Bayes factors for data
tempLogLLH_Ind = zeros(1,n);
tempLogLLH_Dep = zeros(1,n);
for i = 1:n
    if (verbose)
        if (mod(i,25000) == 0)
            fprintf(1, 'Computing Bayes Factors: #%d\n',i);
        end
    end
    if (i == phencol)
        continue;
    end
    if (~disc(i) && disc(phencol))
        % use discrete to continuous likelihood computation
        [tempModel]=learnLocalBN_DiscToCont([], data(:,i)', priorPrecision);    
        tempLogLLH_Ind(i) = tempModel.logLLH;
        [tempModel]=learnLocalBN_DiscToCont(pheno', data(:,i)', priorPrecision);    
        tempLogLLH_Dep(i) = tempModel.logLLH;    
    elseif (disc(phencol) && disc(i))
        % use discrete likelihood computation
        [tempModel]=learnLocalBN_DiscToDisc([], data(:,i)', priorPrecision);    
        tempLogLLH_Ind(i) = tempModel.logLLH;
        [tempModel]=learnLocalBN_DiscToDisc(pheno', data(:,i)', priorPrecision);    
        tempLogLLH_Dep(i) = tempModel.logLLH;
    else
        % use continuous to continuous likelihood computation
        [tempModel]=learnLocalBN_ContToCont([], data(:,i)', priorPrecision);    
        tempLogLLH_Ind(i) = tempModel.logLLH;
        [tempModel]=learnLocalBN_ContToCont(pheno', data(:,i)', priorPrecision);    
        tempLogLLH_Dep(i) = tempModel.logLLH;
    end
end
BFs=tempLogLLH_Dep-tempLogLLH_Ind;




