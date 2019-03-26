function nodes = Triangulate(nodes, eo)
%nodes = Triangulate(nodes, eo)
%
% TRIANGULATE the graph of NODES
% nodes should already be MORALIZED
%
% INPUT : 
%   NODES : list of nodes representing the original bayesian network
%   EO : eliminatino ordering of those nodes
% 
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


for i = 1:length(nodes)
    n = nodes(eo(i));
    % now triangulate n by joining the parents and children -
    % if the link is between nodes yet-to-be-processed according to EO
    % parents should have already been joined in the MORALIZE step
    neighbors = unique([n.pindex, n.cindex],'legacy');
    if (length(neighbors) >= 2)
        for p1=1:length(neighbors)-1
            for p2=p1+1:length(neighbors)
                % only add the edge if the edge is between nodes later in
                % the elimination ordering
                if (find(eo == neighbors(p1)) > i && find(eo == neighbors(p2)) > i)
                    nodes(neighbors(p1)) = addedge(nodes(neighbors(p1)),nodes(neighbors(p2)));
                    nodes(neighbors(p2)) = addedge(nodes(neighbors(p2)),nodes(neighbors(p1)));
                end
            end
        end
    end 
end