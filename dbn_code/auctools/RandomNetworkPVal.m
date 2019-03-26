function [tpval] = RandomNetworkPVal(priorPrecision, data, cols, pheno, numsims, ALG, disc)
%[tpval] = RandomNetworkPVal(priorPrecision, data, cols, pheno, numsims, ALG, disc)
%
% Compares performance of the real CGBN network learned on the data vs.
% performance of the same algorithms on that data with permuted phenotype
% labels.
%
% INPUT:
% PRIORPRECITION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% DATA: data array
% COLS: column names, a cell array of strings
% NUMSIMS: number of permutations to test against.  Uses adaptive
%   permutatino testing and will abandon tests that do not appear to reach
%   promising statistical significance.  Max possible significance is
%   limited by 1/NUMSIMS.
% ALG: Network Learning algorithm indicator (optional)
%   ALG == 1 : use K2 (defualt) search for building networks
%   ALG == 2 : use PhenoCentric search for building networks
%   ALG == 3 : use Exhaustive search for building networks
% DISC: (optional) can specify which columns should be treated as
%   discrete
%
% OUTPUT:
% TPVAL: pvalue.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    ALG = 1;
end
if (nargin < 7)
    disc = IsDiscrete(data);
end

verbose = false;
BFTHRESH = 0;

if (ALG == 3)
    BN = FullBNLearn(data, cols, pheno, BFTHRESH, '', priorPrecision, disc);
    MBNet = BN.MakeIntoMB();
elseif (ALG == 2)
    MBNet = LearnPhenoCentric(data, cols, pheno, priorPrecision, BFTHRESH, verbose, disc);
else
    MBNet = LearnStructure(data, cols, pheno, priorPrecision, '', verbose, disc);
end

if (~isempty(MBNet) && length(MBNet.mb) > 1)
    [realpval] = BNLearnAndTest(MBNet, data, cols);
else
    realpval = 0.5;
end
    
tpval = adaptivepermtester(numsims, realpval, 1, @adaptiverandnetworker, ...
    priorPrecision, data, cols, pheno, ALG, disc);

end
