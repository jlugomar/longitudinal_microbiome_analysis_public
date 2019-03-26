function newfs = factorreduce(cpts, varnames, evidence)
%
% select rows from cpts.table that match VARNAMES on EVIDENCE, wihch are
% corresponding values of the variables identified in VARNAMES.
%
% CPTS: array of class CondProbTable
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


newfs = cpts;

for i = 1:length(cpts)
    newfs(i) = factorreduceworker(cpts(i), varnames, evidence);
end
