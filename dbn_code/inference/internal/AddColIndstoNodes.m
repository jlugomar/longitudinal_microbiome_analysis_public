function nodes = AddColIndstoNodes(nodes, colnames)
% nodes = AddColIndstoNodes(nodes, colnames)
%
% adds or overwrites the exising mapping of nodes to data columns by
% checking each node's NODE.SELF string and matching it with the index of
% NODE.SELF in COLNAMES{}.  This index is stored in the NODE.COLIND field.
%
% Can be used to get a new set of nodes that will test nodes learned from
% one dataset on a different dataset, with a different column order.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

% check that nodes have a self name:
for i = 1:length(nodes)
    if (isempty(nodes(i).self))
        nodes(i).self = colnames(i);
    end
end

% if nodes are given names in node.self, match these to the column names
for i=1:length(colnames)
    for j=1:length(nodes)
        if (strcmpi(colnames(i), nodes(j).self))
            nodes(j).colind = i;
            break;
        end
    end
end
