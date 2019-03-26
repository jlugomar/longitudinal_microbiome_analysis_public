function p = preddistdiffworker(predvalues, datavals)
%
% randomly permute PREDVALUES and then classify DATAVALS based on 
% PREDVALUES > / < 0.5 and then compute the
% p-value of the difference between the two distributions of DATAVALS
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

predvalues = predvalues(randperm(length(predvalues)));

class = predvalues > 0.5;
c1 = datavals(class);
c2 = datavals(~class);
[~, p] = ttest2(c1, c2);

