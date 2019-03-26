function cycle = hasCycle(adjmat)
% cycle = hasCycle(adjmat)
%
% Function checks to see if the adjacency matrix passed in represents a
% directed acyclic graph.
%
% INPUT: 
% ADJMAT: Adjacency matrix representing a bayes net
%
% OUTPUT: 
% CYCLE: true or false, indicating if a cycle was found.
%
% copyright Michael McGeachie 2012.

cycle = false;
for i = 1:length(adjmat)
    cycle = dfsCycleCheck(adjmat, i);
    if (cycle)
        return;
    end
end

    