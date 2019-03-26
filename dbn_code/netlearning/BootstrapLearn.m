function BootsAdjMat = BootstrapLearn(data, cols, pheno, priorPrecision, nboots, ALG, verbose, BFTHRESH, disc)
%BootsAdjMat = BootstrapLearn(data, cols, pheno, priorPrecision, nboots, ALG, verbose, BFTHRESH, disc)
%
% Function to use bootstrapping to generate multiple datasets and learn on
% each a hybrid bayesian network.  
% outputs an adjacency matrix that contains fractional proportions for each
% edge, where the fraction is the proportion of bootstrap realizations that
% contain the edge.
%
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% NBOOTS: number of bootstrap realizations to perform.
% ALG: Network Learning algorithm indicator (optional)
%   ALG == 1 : use K2 (defualt) search for building networks
%   ALG == 2 : use PhenoCentric search for building networks
%   ALG == 3 : use Exhaustive search for building networks
% VERBOSE: optional. if true, will output files of each bootstrap network.
%   Default = false.
% BFTHRESH: optional. Minimum gain in log likelihood required for adding an
%   edge to the network. Default = 0.
% DISC: optional. Boolean for (COOLS). Use to set discrete / continuous
%   variables.
%
% OUTPUT:
% BootsAdjMat: an adjacency matrix which has fractional edges.  Represents the
%   proportion of bootstrap networks which contain that edge.
%
% Copyright Michael McGeachie 2012.  MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    ALG = 1;
end
if (nargin < 7)
    verbose = false;
end
if (nargin < 8)
    BFTHRESH = 0;
end
if (nargin < 9)
    disc = IsDiscrete(data);
end    
adjmat = zeros(length(cols));

% generate a bunch of bootstrap data realizations
for i = 1:nboots
    fprintf(1,'Starting Boot No. %d\n',i);
    rindex = randi(size(data,1),size(data,1),1);
    d = data(rindex,:);
    
    % then call LearnStructure repeatedly
    if (ALG == 2)
        % use phenocentric search with BFTHRESH = 0;
        [~, FullBN] = LearnPhenoCentric(d, cols, pheno, priorPrecision, BFTHRESH, false);
    elseif (ALG == 3)
        FullBN = FullBNLearn(d, cols, pheno, BFTHRESH,  ...
            ['boot',num2str(i),'_exhaustive_'], priorPrecision, disc, verbose);
    else
        % returns a markov blanket BayesNet object
        [~, FullBN] = LearnStructure(d, cols, pheno, priorPrecision, ['boot',num2str(i)], ...
            false);
    end
    % need to reassemble into a adjacency matrix
    adjmat = adjmat + FullBN.adjmat;
end

BootsAdjMat = adjmat ./ nboots;
