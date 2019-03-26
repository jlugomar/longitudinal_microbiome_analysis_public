function tables = GetProbTables(nodes, tree)
%tables = GetProbTables(nodes, tree)
%
% Take a hybrid bayesian network and get the conditional probability tables
% for all the nodes.  Output them in a simple cell array.
%
% INPUT:
% NODES: a NODE array representing the BN
% TREE: the TREE array representing the BN
%
% OUTPUT:
% TABLES: the probability tables, one per discrete node in the BN.
%
% (c) Michael McGeachie 2012.

% have to loop through the nodes, and find each node somewhere in the tree
tables = cell(size(nodes));
for i = 1:length(nodes)
    for j = 1:length(tree)
        % ClusterSetTree objects
        if ((tree(j).index == i) || ...
            (tree(j).discrete && sum (tree(j).dmembers == i) > 0))
            % found the node
            % get the cpt:
            if (tree(j).discrete)
                tables{i} = tree(j).cpt;
            else
                tables{i} = tree(j).lppotential;
            end
            break;
        end
    end
end


