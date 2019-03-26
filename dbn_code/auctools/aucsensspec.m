function [sens, spec] = aucsensspec(confs, truec, cutoff)
% [sens, spec] = aucsensspec(confs, truec, cutoff)
% 
% Computes Sensitivity and Specificity for a particular confidence cuttoff
% in an AUC.
%
% INPUT: 
% CONFS: judgements between [0,1] for membership in a binary class, TRUEC.
% TRUEC: true class for each instance.  Binary class, only.  Parallel array
%   to CONFS.
% CUTOFF: threshold level between [0,1] above which instances will be
%   called class = 1, and below which instances will be called class = 0.
%   Output is based on this CUTOFF value.
%
% OUTPUT:
% SENS: Sensitivity of the given CUTOFF.
% SPEC: Specificity of the given CUTOFF.
%
% 
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


[confs, sorder] = sort(confs);
truec = logical(truec(sorder));


calls = confs >= cutoff;
% false positives:
fp = sum(calls(~truec));
% false negatives
fn = sum(calls(truec) == 0);
% true positives:
tp = sum(calls(truec));
% true negatives:
tn = sum(calls(~truec) == 0);

spec = tn / (tn + fp);

sens = tp / (tp + fn);


