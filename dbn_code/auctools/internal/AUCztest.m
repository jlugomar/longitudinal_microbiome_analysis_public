function [p,z] = AUCztest(diff, var)
%
% this is the simpler test used by the R-code forums post and
% by the Lasko et. al. perl script
%
%

% z statistic
z = diff ./ sqrt(var);

% two sided p value
p = 2 * (1-normcdf(z));
