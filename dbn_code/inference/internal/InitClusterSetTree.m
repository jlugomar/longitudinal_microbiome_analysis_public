function tree = InitClusterSetTree(tree, nodes)
%tree = InitClusterSetTree(tree, nodes)
%
% Initialize the data structures (various indices) of the ClusterSetTree
% passed in the TREE input argument.
%
% Input: 
%   TREE : an array of class ClusterSetTree instantiations.
%   NODES : master array of the BayesNetwork, generated from function
%           NODETREE()
%
% Initailizes several fields in the TREE array:
%   TREE(i).DMEMBERS
%   TREE(i).DVALUES
%   TREE(i).POSTBAG
%   TREE(i).LPPOTENTIAL
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.



%intialize lists
for i = 1:length(tree)
    configs = 1;
    if (~tree(i).discrete)
        % loop through and decide how many members are discrete
        % members, and therefore how many cells are needed in these:
        for j = 1:length(tree(i).members)
            if (nodes(tree(i).members(j)).discrete)
                % could also check the discrete NODE's list of values to
                % calculate the size of the lists
                if (sum(find(tree(i).dmembers == tree(i).members(j))) == 0)
                    tree(i).dmembers(end+1) = tree(i).members(j);
                    % set up dvalues here, too.
                    tree(i).dvalues{end+1} = nodes(tree(i).members(j)).values;
                end
                configs = configs * length(nodes(tree(i).members(j)).values);
            end
        end
    else
        % in a discrete tree, the index itself is part of dmembers:
        tree(i).dmembers = [tree(i).index, tree(i).members];
        tree(i).dvalues = cell(1,length(tree(i).dmembers));
        for j = 1:length(tree(i).dmembers)
            tree(i).dvalues{j} = nodes(tree(i).dmembers(j)).values;
            configs = configs * length(nodes(tree(i).dmembers(j)).values);
        end
    end
    tree(i).lppotential = cell(configs,1);
    tree(i).postbag = cell(configs,1);
end