function adj = AdjMatrixNodes(nodes, disc)
%
% loop through a list of NODES that form a bayesian network and output an
% adjacency matrix.  Will only use the DISCRETE or CONTINUOUS nodes if the
% 2nd argument is specified (true/false), otherwise uses all nodes
%
% INPUT :
%   NODES : array of nodes representing the bayesian network
%   DISC : boolean, if true generates a adj matrix of only the DISCRETE
%       nodes, otherwise only the CONTINUOUS nodes.  Defaults to using ALL
%       nodes when this argument is omitted.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    useall = true;
else
    useall = false;
end

% generate a sparse adjacency matrix for nodes
adj = sparse([],[],[],length(nodes),length(nodes), 3 * length(nodes));
%adj = zeros(length(nodes));
for i = 1:length(nodes)
    if (useall || nodes(i).discrete == disc)
        for j = 1:length(nodes(i).cindex)
            if (useall || nodes(nodes(i).cindex(j)).discrete == disc)
                adj(i,nodes(i).cindex(j)) = 1;
            end
        end
    end
end
end
