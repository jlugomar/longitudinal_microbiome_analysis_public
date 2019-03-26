function [ps, auc] = NBClassifierPvals(data, inds, plotind_vsrandom)
%[ps, auc] = NBClassifierPvals(data, inds, plotind_vsrandom)
%
% Compute the pvalues for difference vs. random guessing with the naive
% bayes classifier based on each (single) attribute indicated by
% data(:,inds) taken one at a time.  Considers each attribute as a
% separate classifier.  Useful for determining which of several attributes
% in a Bayes Net is the best predictor on its own.  Makes comparison to
% random guessing, which is guessing at the baserate of a binary phenotype.
% Makes AUC/ROC plot.
%
% INPUT: 
% DATA: numeric data formatted in columns.  DATA(:,1) must be the PHENOTYPE
%   to predict.
% INDS: columns of DATA to use in a naive-bayes classifier.  Must be
%   discrete variables.
% PLOTIND_VSRANDOM: if provided, the DATA(:,PLOTIND_VSRANDOM) is
%   highlighted on the AUC plot.  If not provided, the max is highlighted.
%
% OUTPUT:
% PS: pvalues vs. random guessing classifier, for each index.
% AUC: AUCs (convex hull) for each index considered as a naive Bayes classifier.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% random guessing according to baserate
baserate = sum(data(:,1)) / length(data(:,1));
t = data(:,1);
base_y = baserate * ones(size(t));
y = ones(length(data(:,1)),length(inds));

ps = ones(size(inds));
auc = ps;
a1 = 0;
for i = 1:length(inds)
    y(:,i) = NaiveBayes2(data,inds(i));
    [ps(i), a1, auc(i), var1, var2, cov] = CompTwoAUC(t, base_y, t, y(:,i), true);
end

if (nargin < 3)
    plotind_vsrandom = find(auc == max(auc),1);
end

figure();
[tp,fp] = roc(t,y(:,plotind_vsrandom));
plot(fp,tp, 'r');
hold on;
[tp,fp] = roc(t,base_y);
plot(fp,tp, 'b');
xlabel('false positive rate');
ylabel('true positive rate');
title('ROC curve');
legend(num2str(auc(plotind_vsrandom)), num2str(a1));
hold off;