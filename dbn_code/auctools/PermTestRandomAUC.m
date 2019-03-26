function [tpval] = PermTestRandomAUC(aucval, data, cols, pheno, numsims)
%[tpval] = PermTestRandomAUC(aucval, data, cols, pheno, numsims)
%
% Compares an AUC on real data vs. random guessing applied to that data.
% Generates random guesses on each DATA point, and computes an AUC.
% Compares real AUC to a population of NUMSIMS random AUCs representing a
% null distribution.
%
% INPUT:
% AUCVAL: the real AUC value, to be tested against random AUCs.
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: string matching one of the COLS names, to use as the phenotype.
% NUMSIMS: number of permutations to test against.  Uses adaptive
%   permutatino testing and will abandon tests that do not appear to reach
%   promising statistical significance.  Max possible significance is
%   limited by 1/NUMSIMS.
%
% OUTPUT:
% TPVAL: pvalue.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

tpval = adaptivepermtester(numsims, aucval, 1, @PermTestRandomAUCWorker,...
    data, cols, pheno);
