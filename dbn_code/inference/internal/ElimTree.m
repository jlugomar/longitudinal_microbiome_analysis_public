function cs = ElimTree(cs, eo)
%cs = ElimTree(cs, eo)
%
% generate an elimination tree with the running interesection property
% from an array of clustersettree objects.
%
% INPUT : 
%   CS : CluserSetTree array of linked objects, a junction tree.
%   EO : elimination ordering for the nodes in the original bayesian
%       network
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% start from the back and move forward:
for i=length(cs):-1:1
    cluster = cs(i);
    high = length(cs);
    parentnode = 0;
    clustereliminated = find(eo == cluster.index);
    % cluster members are nodes, not other clusters:
    for j=1:length(cluster.members)
        melim = find(eo == cluster.members(j));
        if (melim <= high && melim > clustereliminated)
            high = melim;
            % the parent of this cluster should be the cluster
            % corresponding to this node:
            parentnode = cluster.members(j);
        end
    end
    if (parentnode > 0)
        % find cluster corresponding to parentnode:
        for j = i:-1:1
            if (cs(j).index == parentnode)
                cs(i).parent = j;
                cs(j).children{end+1} = i;
                break;
            end
        end
    end    
end

