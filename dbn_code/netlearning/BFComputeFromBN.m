function bfmat = BFComputeFromBN(BN)
% 
% take in a BN of class BayesNet and try to recompute the original Bayes
% Factors that would have given rise to the network.  This is done by
% checking each existing edge, and computing the BF for its inclusion in
% the network based on the source and sink nodes.  Useful for adding edge
% weights to consensus networks obtained from bootstrapping.
%
% INPUT:
% BN: A BayesNet class object.
%
% OUTPUT:
% BFMAT: A matrix matching the adjacency matrix of BN that includes
%   recomputed Bayes Factors for each edge present in BN.adjmat.  Suitable
%   for assignment : BN.weightMatrix = BFMAT;
%
% 
% copyright Michael McGeachie, 2013. MIT license. See cgbayesnets_license.txt.


bfmat = zeros(length(BN.cols));

for i = 1:length(BN.cols)
    
    % for each node
    child = i;
    % find all the parents of this node, these are inds into cols{}
    parents = find(BN.adjmat(:,i) > 0);

    priorparents = [];
    for p = 1:length(parents)
        
        bfs_parents = -inf * ones(1,length(parents));
        for j = 1:length(parents)
            % compute BF for adding this parent to the child for each of the 
            % possible parents
            if (~isempty(intersect(j,priorparents)))
                % if the parent was already added, skip it.
                continue;
            end
            % compute BF
            bfs_parents(j) = bayesfactor_worker(BN.data(:,parents(priorparents)), ...
                BN.data(:,parents(j)), BN.data(:,child), BN.disc(parents(priorparents)), ...
                BN.disc(parents(j)), BN.disc(child), BN.priorPrecision);
        end
        
        % add best parent to prior parents:
        [maxval, bestparentindex] = max(bfs_parents);
        priorparents = [priorparents, bestparentindex];

        % find the parent with the max BF, that would have been added first.
        % take parent with top score, remove it from possible parents.
        bestparent = parents(bestparentindex);
        bfmat(bestparent, child) = maxval;
        
    end
end



