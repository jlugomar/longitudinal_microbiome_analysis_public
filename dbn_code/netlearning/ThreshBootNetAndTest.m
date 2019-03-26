function [convexhullAUC, testconvexhullAUC] = ThreshBootNetAndTest(adjmat, ...
    data, cols, testdata, testcols, pheno, priorPrecision, verbose, thresh)
% convexhullAUC = ThreshBootNetAndTest(adjmat, data, cols, testdata, testcols, pheno, priorPrecision, verbose, thresh)
%
% function to take a continuous adjacency matrix representing the aggregate
% network generated from several bootstrap iterations, using
% BootstrapLearn.m.  Ouputs predictive accuracy, and some trivial graph
% format files of the network.
%
% INPUT:
% ADJMAT: continuous adjacency matrix represent edge frequencies accross
%   bootstrap realizations of a CG Bayes network.
% DATA: training data array
% COLS: column names, a cell array of strings
% TESTDATA: testing data array; must have same columns as DATA
% TESTCOLS: column names, a cell array of strings; must have same columns
%   as data.
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
% testconvexhullAUC: outputs the convex hull AUC of the ROC curve for the
%   aggregate network trained on DATA and tested on TESTDATA.
%
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.


if (nargin < 7)
    thresh = 0.4;
end

disc = IsDiscrete(data);
convexhullAUC = zeros(1,length(thresh));
testconvexhullAUC = zeros(1,length(thresh));

for i = 1:length(thresh)
    t = thresh(i);

    BootstrapNet = GetThreshBootNetwork(adjmat, data, cols, pheno, t, verbose);
    BootstrapNet.priorPrecision = priorPrecision;
    MBNet = BootstrapNet.MakeIntoMB();
    fprintf(1,['Consensus markov-blanket has ', num2str(length(MBNet.mb)), ' nodes\n']);
    MBNet = LearnParamsBN(MBNet);
    
    fprintf(1,'Predicting on Training Data\n');
    [acc, p, z] = PredictPheno(MBNet, verbose);
    convexhullAUC(i) = AUCWorker(acc,p,z,MBNet.GetPhenoCol(),true);

    fprintf(1,'Predicting on Testing Data\n');
    MBNet = MBNet.ReplaceData(testdata,testcols);
    [testacc, testp, testz] = PredictPheno(MBNet, verbose);
    testconvexhullAUC(i) = AUCWorker(testacc,testp,testz,MBNet.GetPhenoCol(),true);

end