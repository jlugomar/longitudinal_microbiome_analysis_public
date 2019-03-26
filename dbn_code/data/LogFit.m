function [auc, good, numvars] = LogFit(data,cols, pheno,verbose)
%[auc, good, numvars] = LogFit(data,cols, pheno,verbose) 
% Fits a logistic regression model, building and testing on data.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% VERBOSE: boolean.  If true, increases output.
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% GOOD: boolean array indictating which variables were included in the
%   model (those with p < 0.05)
% NUMVARS: size of each fold, in number of variables. 
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 4)
    verbose = true;
end

% find pheno col
phencol = strmatch(pheno, cols, 'exact');
xind = true(size(cols));
xind(phencol) = false;
y = data(:,phencol);
x = data(:,xind);


% fit logistic regression w/no 2nd-order interaction terms:
if (verbose)
    fprintf(1,'Learning Logistic Model:\n');
end
[b,dev,stats] = glmfit(x,y,'binomial','link','logit');

% use the values that are statistically significant:
good = stats.p(2:end) <= 0.05;
numvars = sum(good);
if (stats.p(1) <= 0.05)
    const = 'on';
else
    const = 'off';
end

% predict on the data
if (verbose)
    fprintf(1,'Computing Logistic Predictions on Training Data:\n');
end
yzs = glmval(b(stats.p <= 0.05),x(:,good),'logit','constant',const);
if (size(yzs,2) > 1)
    yzs = yzs(:,2);
end
ypred = yzs > .5;
acc = ypred == data(:,phencol);
% compute AUC of these predictions:
auc = AUCWorker(acc,ypred,yzs,data(:,phencol),true,true,verbose);
