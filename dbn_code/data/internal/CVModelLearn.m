function [auc, numnodes, classacc] = CVModelLearn(data, cols, pheno, priorPrecision, ...
    folds, experimentname, verbose)
%[auc, numnodes, classacc] = CVModelLearn(data, cols, pheno, priorPrecision, folds, experimentname, verbose)
% Does crossvalidation of model building and testing on data.  Will build a
% different network for each fold of cross validation.  Provides a
% good test of various parameter settings if called with differen values of
% PRIORPRECISION.
%
% DEPRICATED
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
% FOLDS: Number of folds in the cross-validation to perform.  Default = 5.
% EXPERIMENTNAME: string that will be used in fileoutput names.  Should
%   represent a valid filename
% VERBOSE: boolean.  If true, increases output.
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% NUMNODES: size of each fold, in number of variables. 
% CLASSAC: accuracy per class of the phenotype.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    folds = 5;
end
if (nargin < 6)
    experimentname = 'bayesnet-CV';
end
if (nargin < 7)
    verbose = true;
end


% find pheno col; cound data
phencol = strmatch(pheno, cols, 'exact');
[ncases,ncols] = size(data);

% split data into N-fold CV sets:
r = ceil(rand(1,ncases) * folds);

auc = 0;
numnodes = zeros(1,folds);


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
    for i = 1:length(cols)
        if (d(i))
            discvals{i} = num2cell(unique(data(:,i),'legacy'));
        else 
            discvals{i} = {};
        end
    end

    % learn a BN for this fold of data:
    MBNet = LearnStructure(cvdata, cols, pheno, priorPrecision, [experimentname,'-fold',num2str(k)], verbose);
    numnodes(k) = length(MBNet)-1;
    if (isempty(MBNet))
        % no network learned; no associations worth making
        acc = .5;
        p = ones(size(cvtest,1),1);
        z = p * .5;
    else
        % there was a network
        if (verbose)
            fprintf(1,'Learning Network Parameters\n');
        end
        [tree, nodes] = LearnParams(MBNet, '', cvdata, cols, priorPrecision, discvals);
        if (verbose)
            fprintf(1,'Predicting on Training Data\n');
        end
        [acc, p, z] = PredictPheno(tree, nodes, '', pheno, cvtest, cols, verbose);
    end
    cvPClass = [cvPClass; p];
    cvPredZs = [cvPredZs; z];
    cvTrueClass = [cvTrueClass; cvtest(:,phencol)];
end
    
[auc, classacc] = AUCWorker(acc,cvPClass,cvPredZs,cvTrueClass, true, false, verbose);
    
