function order = TopSort(adjmat)
% simple function to topologically sort a network from the adjacency
% matrix.
% Works on Directed Acyclic Graphs

olist = 1:length(adjmat);
order = zeros(1,length(adjmat));
oind = 1;
added = true;

% this loops until we go once through the network without finding a node to
% add to the topsorted ordering.  This is the case if there's a cycle.
while (added)
    % check each node.
    added = false;
    for i = 1:length(adjmat)
        % if the node is already in the top sort, then skip it.
        if (sum(order == i) == 0)
            parents = olist(adjmat(:,i) > 0);
            % if this has no parents, add it
            if (isempty(parents))
                order(oind) = i;
                oind = oind + 1;
                added = true;
            else
                % otherwise check to see if each parent is in the order
                ok = true;
                for p = 1:length(parents)
                    if (sum(order == parents(p)) == 0)
                        ok = false;
                        break;
                    end
                end
                if (ok)
                    order(oind) = i;
                    oind = oind + 1;
                    added = true;
                end                    
            end
        end
    end
end
