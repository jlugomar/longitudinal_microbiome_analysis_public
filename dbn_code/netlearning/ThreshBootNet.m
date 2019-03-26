function [convexhullAUC, pvp] = ThreshBootNet(adjmat, data, cols, pheno, priorPrecision, verbose, thresh)
% convexhullAUC = ThreshBootNet(adjmat, data, cols, pheno, priorPrecision, verbose, thresh)
%
% function to take a continuous adjacency matrix representing the aggregate
% network generated from several bootstrap iterations, using
% BootstrapLearn.m.  Ouputs predictive accuracy, and some trivial graph
% format files of the network.
%
% INPUT:
% ADJMAT: continuous adjacency matrix represent edge frequencies accross
%   bootstrap realizations of a CG Bayes network.
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
% VERBOSE: if true increases output;
% THRESH: the threshold for including an edge in the aggregate network.
%   Must be in [0,1], used to cut edges from non-edges in the ADJMAT. 
%   Default = 0.4, can be an array of thresholds.
%
% OUTPUT:
% convexhullAUC: outputs the convex hull AUC of the ROC curve for the
%   aggregate network trained and tested on DATA.
% PVP: matrix of p-values for each network computed at a different
%   threshold, vs. each other network computed. 
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.


if (nargin < 7)
    thresh = 0.4;
end

convexhullAUC = zeros(1,length(thresh));

% save phenotype column for AUC calcs:
phncol = find(strcmp(pheno, cols));
phncol = data(:,phncol);

for i = 1:length(thresh)
    t = thresh(i);

    BootstrapNet = GetThreshBootNetwork(adjmat, data, cols, pheno, t, verbose);
    BootstrapNet.priorPrecision = priorPrecision;
    MBNet = BootstrapNet.MakeIntoMB();
    if (verbose)
        MBNet.WriteToGML(['ThreshBootMB_t',num2str(t)]);
    end
    fprintf(1,['Consensus markov-blanket has ', num2str(length(MBNet.nodes)), ' nodes\n']);
    MBNet = LearnParamsBN(MBNet);
    fprintf(1,'Predicting on Training Data\n');
    [acc, p, z] = PredictPheno(MBNet);

    % output AUC or accuracy measure:
    convexhullAUC(i) = AUCWorker(acc,p,z,phncol,true);
    z(p == 0) = 1 - z(p == 0);
    zs(:,i) = z;
end

pvp = -1*ones(length(thresh));
for i = 1:length(thresh)-1
    for j = i+1: length(thresh)
        pvp(i,j) = CompTwoAUC(phncol, zs(:,i), phncol, zs(:,j), true);
    end
end

