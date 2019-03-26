function BNFinal = ForceBuildBN(data, cols, pheno, priorPrecision, onlyMB)
% BNFinal = ForceBuildBN(data, cols, pheno, priorPrecision)
%
% repeatedly call AddNodeToBN() to add each column in COLS to the a Bayes
% Net, regardless of whether or not this increases or decreases the
% posterior likelihood of the data.
%
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
% ONLYMB: optional (default = true), if true, forces each node into the
%   markov blanket of hte phenotype so that each one affects the prediction
%   of the phenotype.
%
% OUTPUT:
% BNFinal: BayesNet Class object representing the obtained bayes net.
%
% copyright Michael McGeachie, 2013.  All Rights Reserved.


% this is a variable to force each node to be added to the MB and have an
% effect on prediction of the phenotype:
if (nargin < 5)
    onlyMB = true;
end
% set up a blank adjmat:
adjmat = zeros(length(cols));
weightMatrix = zeros(length(cols));

% get pheno index so we can skip it:
f = strcmpi(pheno, cols);
phnind = find(f);



for i = 1:length(cols)
    if (i == phnind)
        continue;
    end
    [adjmat, weightMatrix] = AddNodeToBN(adjmat, data, cols, pheno, cols(i), ...
        priorPrecision, onlyMB, weightMatrix);
end

BNFinal = BayesNet([],[],adjmat,weightMatrix,[],true,[],data,cols,pheno,priorPrecision,{});


