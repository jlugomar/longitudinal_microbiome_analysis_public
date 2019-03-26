function [BN, MBnet, mval] = NEdgesFromBoots(adjmat, data, cols, pheno, N, disc)
%[BN, masternetwork, numnodes] = NEdgesFromBoots(adjmat, data, cols, pheno, N)
%
% take a continuous adjacency matrix and order the edges by inclusion
% probability. iteratively try adding each edge until we end up with N
% nodes in the markov blanket of the phenotype.  May have to disallow
% certain edges if they result in cycles.  May result in networks with an
% additional node if adding the final edge results in a Markov Blanket that
% includes two nodes.
%
% INPUT:
% ADJMAT: continuous adjacency matrix representing the edge inclusion
%   frequencies of each edge between each pair of nodes
% DATA: data set matching COLS
% COLS: column names (variable names) for DATA
% PHENO: string of phenotype column in COLS
% N: nubmer of nodes to include in the BN.
%
% OUTPUT:
% BN: BayesNet class representing the network
% MBnet: BayesNet class representing the markov blanket of the network
% MVAL: 
%
% copyright Michael McGeachie, 2013. MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    disc = IsDiscrete(data);
end

masternetwork = zeros(size(adjmat));
weightMatrix = zeros(size(adjmat));
numnodes = 1;
done = false;

% set up return values
BN.adjMatrix = zeros(size(adjmat));
BN.cols = cols;
BN.data = data;


while (~done)
    % Find best edge:
    [~,maxi] = max(max(adjmat,[],2));
    [mval,maxj] = max(adjmat(maxi,:));
    
    if (mval == 0)
        % ran out of edges
        done = true;
        continue;
    end
    
    % add edge to a network
    masternetwork(maxi,maxj) = 1;
    weightMatrix(maxi,maxj) = adjmat(maxi,maxj);

    % and remove this edge from future consideration:
    adjmat(maxi, maxj) = 0;

    
    if (hasCycle(masternetwork))
        % if there was a cycle, drop this edge and continue to the next
        masternetwork(maxi,maxj) = 0;
        continue;
    end
    
    % assign data to a BayesNet structure
    BN = BayesNet([],'mb_bootsnet', masternetwork, weightMatrix, [], false, disc, data, cols, pheno, [],{});
    
    % get markov blanket from the network
    MBnet = BN.MakeIntoMB;

    % number of nodes is minus the phenotype.  So there's an extra node.
    numnodes = length(MBnet.cols);
    
    if (numnodes >= N)
        done = true;
    end
end


