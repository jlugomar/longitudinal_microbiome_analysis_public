function [p, pvalpval] = PredDistDiff(predvalues, datavals, numsims)
%[p, pvalpval] = PredDistDiff(predvalues, datavals, numsims)
%
% Classify datavals based on PREDVALUES > / < 0.5 and then compute the
% p-value of the difference between the two distributions of DATAVALS
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

class = predvalues > 0.5;
c1 = datavals(class);
c2 = datavals(~class);
[~, p] = ttest2(c1, c2);

% now call adaptivepermtester to get a population of these p-values
% when calling adaptivepermtester with a p-value, or something where low is
% better, we just invert the returned value.
pvalpval = adaptivepermtester(numsims, p, 0, @preddistdiffworker, predvalues, datavals);
