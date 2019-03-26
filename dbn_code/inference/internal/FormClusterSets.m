    function clusterset = FormClusterSets(nodes, eo)
% 
% forms a junction tree from the bayes network described by NODES.
%
% INPUT : 
%   NODES : array of nodes (structs with certain properties) arranged into
%       a bayesian network
%   EO : elimination ordering of those nodes
%   
% nodes should be ordered according to the elimination ordering
% however this algorithm proceeds from the end of the array forward.
%
% a CLUSTERSET for a node i is the node i and all neighbor nodes that are
% LATER in the elimination ordering.
%
% OUTPUT : 
%   CLUSTERSET : array of ClusterSetTree class objects forming the junction
%       tree of the NODES
% 
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


clusterset = repmat(ClusterSetTree(), size(eo));

% reverse the elim ordering
eo = eo([end:-1:1]);
for i=length(eo):-1:1
    n = nodes(eo(i));
    clusterset(i).name = n.self;
    clusterset(i).index = n.index;
    neighbors = unique([n.pindex,n.cindex],'legacy');
    % a cluster is discrete if all members are discrete:
    clusterset(i).discrete = n.discrete;
    for j=1:length(neighbors)
        if (find(eo == neighbors(j)) < i)
            clusterset(i).members(end+1) = neighbors(j);
            if (~nodes(neighbors(j)).discrete)
                clusterset(i).discrete = false;
            else
                % keep a separate list of discrete members:
                clusterset(i).dmembers(end+1) = neighbors(j);
                clusterset(i).dvalues{end+1} = nodes(neighbors(j)).values;
            end
        end
    end
end
