function tree = InitClusterPotential(tree, nodes, data, priorPrecision)
%tree = InitClusterPotential(tree, nodes, data, priorPrecision)
%
% Cowell's algorithm 5.1 for "pre-initializing the tree"
%
% Input: 
%   TREE : an array of class ClusterSetTree instantiations.
%   NODES : master array of the BayesNetwork, generated from function
%           NODETREE()
%   DATA : an array of raw data, used to initialize regressions and
%           discrete distributions.
%   COLS : column lables for the DATA input. These must match the
%           NODES(i).SELF strings.
%
% This algorithm primarily populates the TREE structure with all its
% LP-potentials. An LP-POTENTIAL is an instantiation of class LPPotential.
%
% DEBUGGING : could be tested with much larger networks.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% for each node, find a clusterset that contains that node.
for i=1:length(nodes)
    if (nodes(i).discrete)
        % for each discrete node X :
        % find a discrete cluster containing (X union pa(X))
        pax = [nodes(i).index, nodes(i).pindex];
        for j = 1:length(tree)
            if (tree(j).discrete)
                cnodes = [tree(j).index, tree(j).members];
                % checks for subset:
                if (isempty(setdiff(pax,cnodes,'legacy')))
                    % "multiply P(X|pa(X)) into the potential of the cluster"
                    % use discrete factor method:
                    if (isempty(tree(j).cpt.factors))
                        tree(j).cpt = InitDiscreteFactor(tree(j),nodes(i),nodes,data,priorPrecision);
                    else
                        % multiply the current CPT with the new one:
                        newcpt = InitDiscreteFactor(tree(j),nodes(i),nodes,data, priorPrecision);
                        newcpt = joindisjointcpts(newcpt,tree(j).cpt);
                        tree(j).cpt = newcpt;
                    end                       
                    break;
                end
            end
        end
    else
        % continuous vars:
        % find a cluster containing (X union pa(X))
        pax = [nodes(i).index, nodes(i).pindex];
        for j = 1:length(tree)
            cnodes = [tree(j).index, tree(j).members];
            % checks for subset:
            if (isempty(setdiff(pax,cnodes,'legacy')))
                % we have found a cluster containing (X union pa(X))
                tree(j) = InitContFactor(tree(j), nodes(i), nodes, data, priorPrecision);
                break;
            end
        end
    end
end

