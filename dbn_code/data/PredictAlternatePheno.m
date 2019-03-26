function auc = PredictAlternatePheno(BN, newpheno, testdata, testcols)
%auc = PredictAlternatePheno(BN, newpheno, testdata, testcols)
%
% take a BN that was learned on one phenotype and then switch to another
% phenotype to see how well it predicts the other variable.  This is
% something that most users will not want to do.  It is essentially a
% measure of the correlation between the old phenotype and the new one.
%
% INPUT
% BN : BayesNet Class object that was learned from the BN.data to predict
%   BN.pheno
% NEWPHENO : string
% TESTDATA : (optional) alternate testing data set for use after training
%   on DATA.  
% TESTCOLS : (optional) alternate columns matching TESTDATA.
%
% OUTPUT
% AUC: prediction on alternate phenotype by the existing Bayes Net.
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

if (nargin < 3)
    testdata = BN.data;
    testcols = BN.cols;
end


% learn params to go with network
MBNet = LearnParamsBN(BN);

% add the new pheno to the existing dataset by just replacing the old pheno
% data col with the new pheno data col
% operate here on the full data, not hte MB
phenocol = strcmp(testcols, MBNet.pheno);
newphenocol = strcmp(testcols, newpheno);
testdata(:,phenocol) = testdata(:,newphenocol);
MBNet = MBNet.ReplaceData(testdata,testcols);
% important not to change teh phenotype name in cols{}, otherwise it will
% conflict with the node's names

verbose = false;
% test with new phenotype
[acc, p, z] = PredictPheno(MBNet, verbose);

% get auc
[auc, classacc] = AUCWorker(acc,p,z, MBNet.GetPhenoCol(), true, false, verbose);



