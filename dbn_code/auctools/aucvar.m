function [v, varargout] = aucvar(t, y)
% [v, varargout] = aucvar(t, y)
% Computes the variance of the measured AUC for an ROC curve
% uses a method by:
%
% DeLong, E. R.; DeLong, D. M. & Clarke-Pearson, D. L. Comparing the areas
% under two or more correlated receiver operating characteristic curves: a 
% nonparametric approach. Biometrics, Quintiles, Inc., Chapel Hill, North 
% Carolina 27514., 1988, 44, 837-845.
% 
% INPUT:
% t = class of data, 
%       t > 0  :  positive
%       t <= 0 :  negative
%
% y = scores of the classifier, in [0,1]
%
% (c) Michael McGeachie 2008

% process targets
t = t > 0;

% tp = scores for true positive
% tn = scores for true negative
tp = y(t);
tn = y(~t);

pvp = ones(size(tp));
pvn = ones(size(tn));

for i = 1:length(pvp)
    pvp(i) = 1 - (sum(tp(i) > tn) / length(tn));
end

for i = 1: length(pvn)
    pvn(i) = sum(tn(i) > tp) / length(tp);
end

v = var(pvp)/length(tp) + var(pvn)/length(tn);

if (nargout > 1)
    varargout(1) = {pvp};
end
if (nargout > 2)
    varargout(2) = {pvn};
end
