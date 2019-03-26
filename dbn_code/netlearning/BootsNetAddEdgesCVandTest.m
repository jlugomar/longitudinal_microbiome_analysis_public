function [cvAUCs, testingAUCs, numnodes] = BootsNetAddEdgesCVandTest(adjmat, ...
    data, cols, tdata, tcols, priorPrecision, pheno, nodelimit, disc)
%[cvAUCs, numnodes, masternetworks] = BootsNetAddEdgesCVandTest(adjmat, ...
%       data, cols, tdata, tcols, priorPrecision, pheno, nodelimit)
%
% This function explores the results of a bootstrapping-produced continuous
% adjacency matrix.  The adjacency matrix has edge probabilities from a
% bootstrap function such as BootstrapLearn().
%
% Takes the continuous adjacency matrix and sorts and orders the edges.  
% Then iteratively adds each edge, in frequency order, and computes the 
% AUC in Cross-Validation of the resulting network.
%
% INPUT
% ADJMAT: continuous adjacency matrix representing edge frequencies of a
%   bayesian network.  Output from BootstrapLearn(), eg.
% DATA: data set matching COLS
% COLS: column names (variable names) for DATA
% TDATA: test data set matching TCOLS
% TCOLS: test column names (variable names) for TDATA
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% PHENO: string of phenotype column in COLS
% NODELIMIT: max number of nodes to consider.  Default = 20. Optional.
%
% OUTPUT
% CVAUCS: array of AUC values for each network.
% NUMNODES: number of nodes actually included in each network.
% MASTERNETWORKS: cell array of adjacency matrices, for each network.
%
% copyright Michael McGeachie, 2013. MIT license. See cgbayesnets_license.txt.

if (nargin < 8)
    nodelimit = 20;
end
if (nargin < 9)
    disc = IsDiscrete(data);
end
FOLDS = 5;

masternetworks = cell(1,nodelimit);
cvAUCs = zeros(1,nodelimit);
testingAUCs = zeros(1,nodelimit);
numnodes = zeros(1,nodelimit);
for i = 2:nodelimit
    % get a BN with a MB with N total nodes:
    % the BN from this has structure elements for data, cols, disc, and
    % adjmat
    [~, MBnet] = NEdgesFromBoots(adjmat, data, cols, pheno, i, disc);
    MBnet.priorPrecision = priorPrecision;
    numnodes(i) = length(MBnet.mb);
    masternetworks{i} = MBnet;

    % learn this in CV 
    [cvAUCs(i), ~] = CVParamTest(MBnet, FOLDS, ['BN_EdgeCycle', num2str(i)], false);
    % test it in test data:
    testingAUCs(i) = BNLearnAndTest(MBnet, tdata,tcols);
end


