function auc = adaptivelogisticpredworker(x, y)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

[nrows,nsnps] = size(x);


r = randperm(nrows);
y = y(r);

warning off;
[b,~,~] = glmfit(x,y,'binomial','link','logit');
warning on;

% predict on the data
yzs = glmval(b,x,'logit');
if (size(yzs,2) > 1)
    yzs = yzs(:,2);
end
ypred = yzs > .5;
acc = ypred == y;
% compute AUC of these predictions:
auc = AUCWorker(acc,ypred,yzs,y,true,true,false);

