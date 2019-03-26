function [auc, classacc] = CVParamTest(MBNet, folds, experimentname, ...
    verbose, randfolds)
%[auc, classacc] = CVParamTest(MBNet, folds, experimentname, verbose, randfolds)
%
% Does crossvalidation of a given network, testing it on data.  Provides a
% good test of the sensitivity of the network's parameters - marginal 
% distributions at each node - to outliers in the data.  Does cross
% validation only on the LearnParams() function.
%
% INPUT:
% MBNET: BayesNet class object representing the markov blanket of the
%   bayesian network. Must include data, cols, pheno, disc, priorPrecision,
%   etc.
% EXPERIMENTNAME: string that will be used in fileoutput names.  Should
%   represent a valid filename
% VERBOSE: boolean.  If true, increases output.
% RANDFOLDS: boolean, if true, will randomly select each cross validation
%   fold.  not gauranteed to be the same size.  (default = true)
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% CLASSAC: accuracy per class of the phenotype.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    folds = 5;
end
if (nargin < 3)
    experimentname = [MBNet.title, '-CV'];
end
if (nargin < 4)
    verbose = true;
end
if (nargin < 5)
    randfolds = true;
end

% find pheno col; count data
data = MBNet.data;
[ncases,ncols] = size(data);

if (randfolds)
    % randomly split data into N-fold CV sets:
    r = ceil(rand(1,ncases) * folds);
else
    r = zeros(1,ncases);
    for i = 1:folds
        r(ceil(ncases/folds)*(i-1)+1:ceil(ncases/folds)*i) = i;
    end
    r = r(1:ncases);
end
auc = 0;


cvPClass = [];
cvPredZs = [];
cvTrueClass = [];
% for each fold in teh cross-validation
for k = 1:folds
    cvdata = data(r ~= k,:);
    cvtest = data(r == k,:);
    if (verbose)
        fprintf(1,'Starting Fold %d!\n',k);
    end
    if (isempty(cvtest))
        if (verbose)
            fprintf(1,'Skipping Fold %d because there is no test data!\n',k);
        end
        continue;
    end
    
    % check values of discrete vars::
    d = IsDiscrete(data, 5);
    discvals = cell(size(d));
    for i = 1:ncols
        if (d(i))
            discvals{i} = num2cell(unique(data(:,i),'legacy'));
        else 
            discvals{i} = {};
        end
    end
    MBNet = MBNet.AddDiscVals(discvals);
    
    % learn params this fold of data:  
    MBNet.data = cvdata;
    if (verbose)
        fprintf(1,'Learning Network Parameters\n');
    end
    MBTestNet = LearnParamsBN(MBNet);
    if (verbose)
        fprintf(1,'Predicting on Training Data\n');
    end
    MBTestNet = MBTestNet;
    MBTestNet.data = cvtest;
%   MBTestNet = MBNet.ReplaceData(cvdata,cols);
    [acc, p, z] = PredictPheno(MBTestNet, verbose);
    cvPClass = [cvPClass; p];
    cvPredZs = [cvPredZs; z];
    cvTrueClass = [cvTrueClass; MBTestNet.GetPhenoCol()];
end
    
[auc, classacc] = AUCWorker(acc,cvPClass,cvPredZs,cvTrueClass, true, false, verbose);
    
