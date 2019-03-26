function [auc, pval, tauc, tpval] = ConditionalLogisticRegression(...
    x, y, conds, tx, ty, tconds, numsims)
%[auc, pval, tauc, tpval] = ConditionalLogisticRegression(x, y, conds, tx, ty, tconds, numsims)
%
% INPUT:
%  X : matrix of independent variables in separate columns
%  Y : binary column vector of class labels.  Must be (0,1)
%  CONDS : matrix of discrete conditioning variables for building separate
%   classifiers according to each value of CONDS
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 7)
    numsims = 100;
end

% get the base AUC on the fitted data:
auc = adaptivecondlogregworker(x,y,conds, false);

% and compute a p-value for this AUC using permutation testing
pval = adaptivepermtester(numsims, auc, 1, @adaptivecondlogregworker, ...
    x, y, conds, true);

if (nargin >= 6)
    % see how the original model does on some test data:
    tauc = adaptivecondlogregworker(x,y,conds,false,tx,ty,tconds);

    % and compute a p-value for this AUC using permutation testing
    tpval = adaptivepermtester(numsims, auc, 1, @adaptivecondlogregworker, ...
        x, y, conds, true, tx, ty, tconds);
else
    tauc = -1;
    tpval = -1;
end


