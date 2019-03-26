function [auc, tpval, npreds] = TrainAndTestNB(traindata,testdata,cols,pheno, BFTHRESH, priorPrecision, varlim)
% [auc, tpval, npreds] = TrainAndTestNB(traindata,testdata,cols,pheno, BFTHRESH, priorPrecision)
%
% Wrapper function to call NBBayesNet(), get a naive-bayes bayesian
% network, and then:
% 1) write to GML
% 2) learn parameters using the TRAINDATA only
% 3) predict on TESTDATA only
% 4) report AUC on TESTDATA
% 5) use PermTestRandomAUC() to report a p-value of that AUC's difference
%       from random prediction
%
%
% INPUT:
%   TRAINDATA: data array (data points in rows by variables in columns) for
%       learning the naive-bayes model 
%   TESTDATA: data array for testing the naive-bayes model. must have same
%       columns in same order as TRAINDATA
%   COLS : column names, a cell array of strings
%   PHENO : a string representing the phenotype column to predict.  Is matched
%     against the COLS array.
%   BFTHRESH : minimum threshold (in log Bayes Factor) for including a
%   PRIORPRECISION : a structure including the usual HybridBayesNets
%     parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%   VARLIM : an upper bound on the number of variables included in the
%       model.  If sum(bayesfactor(TRAINDATA(:,k)) > BFTHRESH) > VARLIM, on
%       the top VARLIM columns will be included in the BayesNet.
%
% OUTPUT : 
%   AUC : AUC of naive-bayes model on TESTDATA
%   TPVAL : permutation p-value for AUC difference from random performance
%   NPREDS : the number of predictive variables found and used in the
%       Naive-Bayes model
%
% Copyright Michael McGeachie, 2015.  MIT license. See cgbayesnets_license.txt.

if (nargin < 7)
    varlim = inf;
end

verbose = false;
numsims = 1000;


fprintf(1,'Learning (NB) naive bayes for predicting %s\n',pheno);
% returns a markov-blanket reduced network:
disc = IsDiscrete([traindata;testdata]);
NBNet = NBBayesNet(double(traindata), cols, pheno, BFTHRESH, ['naivebayes_',pheno], ...
    priorPrecision, disc, verbose, varlim);
%NBNet.WriteToGML();
npreds = length(NBNet.mb) - 1;
if (npreds == 0)
    fprintf(1,'No predictors found!\n');
    auc = .5;
    tpval = 1;
    return;
end
% test NB
fprintf(1,'Learning NB network distribution parameters\n');
NBNet = NBNet.LearnParams();
%
fprintf(1,'Predicting NB on Test Data for %s\n',pheno);
NBNettest = NBNet.ReplaceData(double(testdata),cols);
[acc, p, z] = NBNettest.Predict(verbose);
[auc, classacc] = AUCWorker(acc,p,z,NBNettest.GetPhenoCol());
[tpval] = PermTestRandomAUC(auc, double(testdata), cols, pheno, numsims);


