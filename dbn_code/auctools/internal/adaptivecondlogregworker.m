function auc = adaptivecondlogregworker(x,y,conds,randomize,tx,ty,tconds)
%
% INPUT:
%  X : matrix of independent variables in separate columns
%  Y : binary column vector of class labels.  Must be (0,1)
%  CONDS : matrix of discrete conditioning variables for building separate
%  classifiers according to each value of CONDS
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 4)
    randomize = true;
end
if (nargin < 5)
    test = false;
else
    test = true;
end

[nrows, condvars] = size(conds);

if (randomize)
    r = randperm(nrows);
    y = y(r);
    if (test)
        r = randperm(length(ty));
        ty = ty(r);
    end
end

% parse conditions:
cvals = cell(1,condvars);
for i = 1:condvars;
    cvals{i} = unique(conds(:,i),'legacy');
end

if (~test)
    acc = -1 * ones(nrows,1);
else
    acc = -1 * ones(length(ty));
end
ypred = acc;
yscores = acc;
% figure out how many different regressions we have:
[cpattern, vinds, numvals] = ValueIncrement(cvals, []);
for i = 1:numvals
    % find matching rows:
    creppat = repmat(cpattern, nrows, 1);
    matches = conds - creppat;
    % matches are those rows who's abs(sum) is zero
    zs = sum(abs(matches),2);
    xsubset = x(zs == 0);
    ysubset = y(zs == 0,:);
    
    % now do a linear model on this subset:
    warning off;
    [b,~,~] = glmfit(xsubset,ysubset,'binomial','link','logit');
    warning on;

    if (test)
        creppat = repmat(cpattern, length(ty),1);
        matches = tconds - creppat;
        % overwrite data subsets here because we only need the test data
        % for predictions in the following
        zs = sum(abs(matches),2);
        xsubset = tx(zs == 0);
        ysubset = ty(zs == 0);
    end

    % get predictions from this model:
    yzs = glmval(b,xsubset,'logit');
    if (size(yzs,2) > 1)
        yzs = yzs(:,2);
    end
    yscores(zs == 0) = yzs;
    ypred(zs == 0) = yscores(zs == 0) > .5;
    acc(zs == 0) = ypred(zs == 0) == ysubset;

    [cpattern, vinds] = ValueIncrement(cvals, cpattern, vinds); 
end

%compute AUC from all predictions
if (~test)
    auc = AUCWorker(acc,ypred,yscores,y,true,true,false);
else
    auc = AUCWorker(acc,ypred,yscores,ty,true,true,false);
end





