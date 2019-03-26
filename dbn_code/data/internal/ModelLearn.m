function [auc, numnodes] = ModelLearn(data, cols, pheno, priorPrecision, ...
   experimentname, verbose)
%[auc, numnodes] = ModelLearn(data, cols, pheno, priorPrecision, experimentname, verbose)
% Build a bayesian network model of the data and then measure prediction on
% that data.  
%
% Depricated.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% EXPERIMENTNAME: string that will be used in fileoutput names.  Should
%   represent a valid filename
% VERBOSE: boolean.  If true, increases output.
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% NUMNODES: size of each fold, in number of variables. 
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    verbose = true;
end
if (nargin < 5)
    experimentname = 'bayesnet-fittedvalue';
end

% find pheno col
phencol = strmatch(pheno, cols, 'exact');

% learn a BN for this fold of data::
MBNet = LearnStructure(data, cols, pheno, priorPrecision, experimentname, verbose);
numnodes = length(MBNet.mb)-1;
if (isempty(MBNet))
    % no network learned; no associations worth making
    acc = .5;
    p = ones(size(data,1),1);
    z = p * .5;
else
    % there was a network
    if (verbose)
        fprintf(1,'Learning Network Parameters\n');
    end
    [tree, nodes] = LearnParams(MBNet, '', data, cols, priorPrecision);
    if (verbose)
        fprintf(1,'Prediction on Training Data:\n');
    end
    [acc, p, z] = PredictPheno(tree, nodes, '', pheno, data, cols, false);
end

auc = AUCWorker(acc,p,z,data(:,phencol));
    
