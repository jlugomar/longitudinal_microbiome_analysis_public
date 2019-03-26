function [auc, numnodes, testauc, MBNet, classacc, testclassacc] = ...
    ModelLearnAndTest(data, cols, testdata, testcols, pheno, priorPrecision, ...
        experimentname, verbose, discnames, ALG)
%[auc, numnodes, testauc, model, classacc, testclassacc] = ModelLearnAndTest(data, cols, testdata, testcols, pheno, priorPrecision, experimentname, verbose, discnames, ALG)
%
% Build a bayesian network model of the data and then measure prediction on
% a test data set, separate from the original training dataset.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% TESTDATA: testing data array
% TESTCOLS: testing column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% PRIORPRECISION: a structure including the usual CGBayesNets parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% EXPERIMENTNAME: string that will be used in fileoutput names.  Should
%   represent a valid filename
% VERBOSE: boolean.  If true, increases output.
% DISCNAMES: cell string array of names of columns to be considered
%   discrete.  These will have extreme values truncated when they do not
%   occur in both the testing and training datasets.
% ALG: Network Learning algorithm indicator (optional)
%   ALG == 1 : use K2 (defualt) search for building networks
%   ALG == 2 : use PhenoCentric search for building networks
%   ALG == 3 : use Exhaustive search for building networks
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% NUMNODES: size of each fold, in number of variables. 
% TESTAUC: AUC of model on test dataset.
% MBNET: Class BayesNet object representing the markov blanket of the the 
%   Bayes Net that was learned.
% CLASSACC: accuracy of model on training data, per class of the phenotype.
% TESTCLASSACC: accuracy of model on test data, per class of the phenotype.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.
%

if (nargin < 8)
    verbose = true;
end
if (nargin < 7)
    experimentname = 'bayesnet-fittedvalue';
end
if (nargin < 9 || isempty(discnames))
    discnames = cols(IsDiscrete(data));
end
if (nargin < 10)
    ALG = 1;
end

DOFILTERING = false;

% find pheno col
phencol = strmatch(pheno, cols, 'exact');
testphencol = strmatch(pheno, testcols, 'exact');

% find other discrete columns:
disc_train = false(1,length(cols));
disc_test = false(1,length(cols));
discvals_train = cell(1,length(cols));
discvals_test = cell(1,length(cols));
for i = 1:length(discnames)
    ind = strmatch(discnames(i), cols, 'exact');
    if (IsDiscrete(data(:,ind)))
        testind = strmatch(discnames(i), testcols,'exact');
        if (~isempty(testind))
            trainvals = unique(data(:,ind),'legacy');
            testvals = unique(testdata(:,testind),'legacy');
            disc_train(ind) = true;
            disc_test(testind) = true;
            if (DOFILTERING)            
                rares = setdiff(testvals,trainvals,'legacy');
                if (~isempty(rares))
                    testdata(testdata(:,testind) == rares(1),testind) = max(trainvals);
                end
            else
                discvals_train{ind} = num2cell(union(testvals, trainvals, 'legacy'));
                discvals_test{testind} = num2cell(union(testvals, trainvals, 'legacy'));
            end
        end
    end
end



% learn a BN for the training data::
BFTHRESH = 0;
if (ALG == 3)
    BN = FullBNLearn(data, cols, pheno, BFTHRESH, experimentname, priorPrecision, disc_train);
    MBNet = BN.MakeIntoMB();
elseif (ALG == 2)
    MBNet = LearnPhenoCentric(data, cols, pheno, priorPrecision, BFTHRESH, verbose, disc_train);
else
    MBNet = LearnStructure(data, cols, pheno, priorPrecision, experimentname, verbose, disc_train);
end

MBNet = MBNet.AddDiscVals(discvals_train);
numnodes = length(MBNet.mb)-1;
if numnodes > 1
    if (verbose)
        fprintf(1,'\t Identified %d variables for predicting %s.\n', numnodes, pheno);
    end
else
    if (verbose)
        fprintf(1,'No variable identified for prediction. Run again after adjusting parameters.\n');
    end
    auc = .5;
    numnodes = 0;
    testauc = .5;
    MBNet = [];
    classacc = [];
    testclassacc = [];     
    return;
end

if (numnodes <= 1)
    % no network learned; no associations worth making
    acc = .5;
    p = ones(size(data,1),1);
    z = p * .5;
    testacc = acc;
    testp = ones(size(testdata,1),1);
    testz = testp * .5;
else
    % there was a network
    if (verbose)
        fprintf(1,'Learning Network Parameters\n');
    end
    % include disc values here:
    MBNet = LearnParamsBN(MBNet);
    if (verbose)
        fprintf(1,'Prediction on Training Data:\n');
    end
    [acc, p, z] = PredictPheno(MBNet);
    %testnodes = AddColIndstoNodes(nodes, testcols);
    MBTestNet = MBNet.ReplaceData(testdata,testcols);
    if (verbose)
        fprintf(1,'Prediction on Testing Data:\n');
    end
    [testacc, testp, testz] = PredictPheno(MBTestNet);
    if (verbose)
        MBNet.WriteToGML();
    end
end
[auc, classacc] = AUCWorker(acc,p,z,MBNet.GetPhenoCol());
[testauc, testclassacc] = AUCWorker(testacc,testp,testz,MBTestNet.GetPhenoCol());
    
