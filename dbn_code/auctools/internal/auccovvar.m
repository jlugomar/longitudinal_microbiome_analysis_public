function cv = auccovvar(t1, y1, t2, y2)
% cv = auccovvar(t1, y1, t2, y2)
% computes covariance of two ROC curves
%
% (c) Michael McGeachie 2008

% get variances
% and "placement values" for positive and negative examples
[v1,pv1,nv1] = aucvar(t1,y1);
[v2,pv2,nv2] = aucvar(t2,y2);

s1 = pv1' * pv2 - length(pv1) * mean(pv1) * mean(pv2);

s2 = nv1' * nv2 - length(nv1) * mean(nv1) * mean(nv2);

cv = s1 / (length(pv1) * (length(pv1) -1)) + s2 / (length(nv1) * (length(nv1) -1));
