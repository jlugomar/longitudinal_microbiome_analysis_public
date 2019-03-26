function [auc, pval, acc, tauc, tpval, tacc, pseudo_r2, dev] = logisticpredworker(x, y, tx, ty, numsims)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

warning off;
[b,dev,~] = glmfit(x,y,'binomial','link','logit');
[~,dev_base,~] = glmfit(ones(size(y)),y,'binomial','link','logit','constant','off');
warning on;

% compute McFadden's pseudo R^2
pseudo_r2 = 1 - dev / dev_base;

[auc, pval, acc] = lpworkerhelper(b, x, y, numsims);
[tauc, tpval, tacc] = lpworkerhelper(b, tx, ty, numsims);

end

function [auc, pval, acc] = lpworkerhelper(b, x, y, numsims)
% predict on the data
yzs = glmval(b,x,'logit');
if (size(yzs,2) > 1)
    yzs = yzs(:,2);
end
ypred = yzs > .5;
acc = ypred == y;
acc = sum(acc)/length(y);
% compute AUC of these predictions:
[auc,acc] = AUCWorker(acc,ypred,yzs,y,true,true,false);

% and compute a p-value for this AUC using permutation testing
pval = adaptivepermtester(numsims, auc, 1, @adaptivelogisticpredworker, ...
    x, y);
end

