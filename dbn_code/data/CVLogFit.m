function [auc, good, numvars, classacc] = CVLogFit(data, cols, pheno, folds, verbose)
%[auc, good, numvars, classacc] = CVLogFit(data, cols, pheno, folds, verbose)
% Does crossvalidation of logistic regression model building and testing on data.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% FOLDS: Number of folds in the cross-validation to perform.  Default = 5.
% VERBOSE: boolean.  If true, increases output.
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% GOOD: boolean array indictating which variables were included in the
%   model (those with p < 0.05)
% NUMVARS: size of each fold, in number of variables. 
% CLASSAC: accuracy per class of the phenotype.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 4)
    folds = 5;
end
if (nargin < 5)
    verbose = true;
end

% find pheno col
phencol = strmatch(pheno, cols, 'exact');
xind = true(size(cols));
xind(phencol) = false;


% split data into N-fold CV sets:
[ncases,ncols] = size(data);
r = ceil(rand(1,ncases) * folds);

auc = 0;
numvars = zeros(1,folds);

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

    % fit logistic regression w/no 2nd-order interaction terms:
    if (verbose)
        fprintf(1,'Learning Logistic Model:\n');
    end
    [b,dev,stats] = glmfit(cvdata(:,xind),cvdata(:,phencol),'binomial','link','logit');

    % use the values that are statistically significant:
    good = stats.p(2:end) <= 0.05;
    numvars(k) = sum(good);
    if (stats.p(1) <= 0.05)
        const = 'on';
    else
        const = 'off';
    end

    % predict on the data
    if (verbose)
        fprintf(1,'Computing Logistic Predictions on Training Data:\n');
    end
    xtest = cvtest(:,xind);
    yzs = glmval(b(stats.p <= 0.05),xtest(:,good),'logit','constant',const);
    if (size(yzs,2) > 1)
        yzs = yzs(:,2);
    end
    ypred = yzs > .5;
    acc = ypred == cvtest(:,phencol);
    % compute AUC of these predictions:
    
    cvPClass = [cvPClass; ypred];
    cvPredZs = [cvPredZs; yzs];
    cvTrueClass = [cvTrueClass; cvtest(:,phencol)];
end

[auc,classacc] = AUCWorker(acc,cvPClass,cvPredZs,cvTrueClass,true,true, verbose);
    
