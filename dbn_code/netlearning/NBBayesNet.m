function outnet = NBBayesNet(data, cols, pheno, BFTHRESH, ...
    outfilename, priorPrecision, disc, verbose, varlim)
% outnet = NBBayesNet(data, cols, pheno, BFTHRESH, outfilename, priorPrecision, disc, verbose, varlim)
% 
% Builds a naive bayes model for predicting the phenotype from the data;
% use a BayesNet class object to compute prediction results.  Uses sparse()
% matrix representation for the adjaceny matrix and weight matrix in the
% BayesNet class object; this represents a slight modification of the
% typical BN objects.
%
% Should be able to deal with a few thousand nodes in the network (ie, with
% bayesfactor > BFTHRESH), out of a total of many thousands of possible
% variables.
%
% INPUT:
%   DATA : data array (data points in rows by variables in columns)
%   COLS : column names, a cell array of strings
%   PHENO : a string representing the phenotype column to predict.  Is matched
%     against the COLS array.
%   BFTHRESH : minimum threshold (in log Bayes Factor) for including a
%     column in the NaiveBayes model.
%   OUTFILENAME : filename for printing out the network structure in a .tgf
%     file.  (Trivial Graph Format).
%   PRIORPRECISION : a structure including the usual HybridBayesNets
%     parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%   DISC : (optional) logical array parallel to cols specifying which are
%       discrete.
%   VERBOSE : increases output printed to stdout. (optional, default =
%     false).
%   VARLIM : an upper bound on the number of variables included in the
%       model.  If sum(bayesfactor(TRAINDATA(:,k)) > BFTHRESH) > VARLIM, on
%       the top VARLIM columns will be included in the BayesNet.
%
% OUTPUT : 
%   OUTNET : BayesNet class object that embodies the Naive Bayes model 
%       identified.
%
%
% Copyright Michael McGeachie, 2015.  MIT license. See cgbayesnets_license.txt.


if (nargin < 7)
    disc = IsDiscrete(data);
end
if (nargin < 8)
    verbose = false;
end
if (nargin < 9)
    varlim = inf;
end

% compute BayesFactors
bfs = BayesFactorScore(data,cols,pheno,priorPrecision,disc,verbose);
netbfs = bfs > BFTHRESH;

% if there's more variables here passing the threshold, just take the top
% VARLIM of them:
if (sum(netbfs) > varlim)
    [~,sorder] = sort(bfs,'descend');
    netbfs = false(size(bfs));
    netbfs(sorder(1:varlim)) = true;
end

% take all cols with BF > BFTHRESH and include in a model
% string match phenotype vs. columns:
match = strcmp(pheno,cols);
phncol = find(match);


inds = 1:length(bfs);
% allocates a sparse adjacency matrix that has links from pheno to
% everything with BF > 0, and also has enough space for a full NB model.
adjmat = sparse(phncol * ones(1,sum(netbfs)), inds(netbfs),ones(1,sum(netbfs)),length(inds),length(inds),length(inds));
weightmatrix = sparse(phncol * ones(1,sum(netbfs)), inds(netbfs),bfs(netbfs),length(inds),length(inds),length(inds));
% make the BayesNet class object
nbbn = BayesNet([], ... % BN object has no node tree yet
    ['nb',outfilename], ... % title of the network
    adjmat, ... % sparse adjacency matrix
    weightmatrix, ... % sparse weight matrix is just the Bayes Factors 
    inds(netbfs), ... % the Markov Blanket is the whole network
    false, ... % it isn't a Markov Blanket network until its reduced to just the included columns
    disc, data, cols, pheno, priorPrecision); 

% convert to Markov Blanket:
outnet = nbbn.MakeIntoMB();

% done





