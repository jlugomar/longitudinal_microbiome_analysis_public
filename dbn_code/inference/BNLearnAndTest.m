function [auc, p, z] = BNLearnAndTest(BN, testdata, testcols)
% [auc, p, z] = BNLearnAndTest(BN, testdata,testcols)
%
% Takes a BayesNet class object that represents a bayesnet with a
% determined structure.  This function learns the parameters for it, using
% the data contained in the BN object, and then tests the BN on predicting
% the phenotype.  Optionally tests on different data than parameters were 
% learned with.
%
% INPUT:
% BN: BayesNet class object.  Must have structure defined and be a
%   markov-blanket network.
% TESTDATA: testing data array (optional).  If empty, will test on 
%   BN.data
% TESTCOLS: testing column names, a cell array of strings (optional). 
%
% OUTPUT:
% AUC: the predictive accuracy of BN on the data
% P: is predicted class value for each testing sample
% Z: is the predicted value (on [0,1]) for each testing sample
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


if (nargin < 2)
    SEPARATE_TEST = false;
else
    SEPARATE_TEST = true;
end


BN = LearnParamsBN(BN);
if (SEPARATE_TEST)
    BNTest = BN.ReplaceData(testdata,testcols);
    [acc, p, z] = PredictPheno(BNTest);
    [auc, classacc] = AUCWorker(acc,p,z,BNTest.GetPhenoCol());
else
    [acc, p, z] = PredictPheno(BN);
    [auc, classacc] = AUCWorker(acc,p,z,BN.GetPhenoCol());
end
