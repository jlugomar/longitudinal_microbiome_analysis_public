function s = logplus(v)
% sum a vector V of logs.
% This is useful when summing a group of log-probabilities
% uses an exact computation.
%
% Sum(v) = m + log(sum(exp(log(v) - m))),
% where m = max (log (v))
%
% 
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

m = max(v);

s1 = sum(exp(v - m));
if (s1 <= 0)
    error('error in LOGPLUS! log of negative values');
end
s = m + log(s1);