function adjmat = AdjMatfromBN(nodes, nodenames)
%adjmat = AdjMatfromBN(nodes, nodenames)
%
% Computes an adjacency matrix from a bayesian network structure.
%
% INPUT: 
% NODES: an array of NODE structures, 
% NODENAMES: cell array of strings of names of those nodes.
%
% OUTPUT: 
% ADJMAT: is ordered according to the column/node names listed in NODENAMES.
%
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.

adjmat = zeros(length(nodenames));

namemap = zeros(1,length(nodenames));
for i = 1:length(nodes)
    hits = strcmpi(nodes(i).self, nodenames);
    namemap(i) = find(hits > 0);
end

for i = 1:length(nodes)
    if (~isempty(nodes(i).children))
        for j = 1:length(nodes(i).children)
            adjmat(namemap(i), namemap(nodes(i).cindex(j))) = 1;
        end
    end
end

