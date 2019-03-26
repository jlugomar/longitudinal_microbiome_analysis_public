function [p, a1, a2, var1, var2, cov] = CompTwoAUC(t1, y1, t2, y2, useCH)
%[p, a1, a2, var1, var2, cov] = CompTwoAUC(t1, y1, t2, y2, useCH)
%
% Computes comparisons between two AUCs computed on the same data.
% Provides a p-value for the difference between them.  Interestingly, the
% data should be the same in t1 and t2, but they do not have to be in the
% same order.
%
% INPUTS:
% y1: is predicted class value for ROC1
% t1: is true class value for ROC1
% y2: is predicted class value for ROC2
% t2: is true class value for ROC2. t1 should == t2
% useCH: if true will do comparisons based on the convex hull of the ROCs
%   involved. (default = false)
%
% OUTPUTS:
% p = p-value for difference between ROC1 and ROC2. Lower p implies greater
%   statistical evidence for difference.
% a1 = AUC of the ROC defined by (y1,t1).
% a2 = AUC of the ROC defined by (y2,t2).
% var1 = estimated variance of AUC1, computed using DeLong et al.'s method.
% var2 = estimated variance of AUC2.
% cov = covariance of the two AUCs
%
% (c) Michael McGeachie 2008

if (nargin < 5)
    useCH = false;
end

% get AUC for ROC curves:
if (useCH)
    % use the convex hull for theoretical opportunism
    [tp1,fp1] = rocch(t1,y1);
    [tp2,fp2] = rocch(t2,y2);
else
    [tp1,fp1] = roc(t1,y1);
    [tp2,fp2] = roc(t2,y2);
end
a1 = auroc(tp1,fp1);
a2 = auroc(tp2,fp2);


% variance of AUCs
var1 = aucvar(t1,y1);
var2 = aucvar(t2,y2);

% compute covariance of AUC
cov = auccovvar(t1,y1,t2,y2);

% variance of the difference
vardiff = var1 + var2 - 2 * cov;

% acutal difference
diff = abs(a1 - a2);

% and a pvalue computed using a z-statistic, assuming a normal distribution
% z statistic
z = diff ./ sqrt(vardiff);
% two sided p value
p = 2 * (1-normcdf(z));
