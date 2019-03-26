function nodes = Moralize(nodes)
%nodes = Moralize(nodes)
%
% take a DAG graph given as an array of NODES and join the parents of every
% node with an edge.  edges are added by adding a link in NODE.CHILDREN
%
% INPUT : 
%   nodes : array of class NODE() representing a bayesian network
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


for i=1:length(nodes)
    if (length(nodes(i).pindex) >= 2)
        for p1=1:length(nodes(i).pindex)-1
            for p2=p1+1:length(nodes(i).pindex)
                nodes(nodes(i).pindex(p1)) = addedge(nodes(nodes(i).pindex(p1)),...
                    nodes(nodes(i).pindex(p2)));
                nodes(nodes(i).pindex(p2)) = addedge(nodes(nodes(i).pindex(p2)),...
                    nodes(nodes(i).pindex(p1)));
            end
        end
    end
end

