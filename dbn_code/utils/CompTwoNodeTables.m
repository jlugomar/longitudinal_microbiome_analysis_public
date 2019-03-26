function same = CompTwoNodeTables(tables1, tables2)
%same = CompTwoNodeTables(tables1, tables2)
%
% compare the Tables output by two calls to different GetProbTables()
% both tables must have the same nodes in the same order. 
% this is useful for determining if the nodes all point in the same
% direction, when the parameters are learned on somewhat different data.
%
% BETA
%
%
% copyright Michael McGeachie 2013.  MIT license. See cgbayesnets_license.txt.

same = true(size(tables1));
% step through each node
% node 1 should be the head node... but possibly some other node

for i = 1:length(tables1)
    nodeclass = class(tables1{i});
    if (strcmpi(nodeclass, 'CondProbTable'))
        % discrete node:
        continue;
    end
    % if the node is a continuous node with parents:
    t1 = tables1{i};
    t2 = tables2{i};
    % should be a PHNVALSx1 cell of LPPotentials
    for j = 1:length(t1)
        lp1 = t1{j};
        lp2 = t1{j};
        % check that params have the same sign
        if (~isempty(lp1.params))
            for k = 1:length(lp1.params)
                same(i) = same(i) && ((lp1.params(k) >= 0) == (lp2.params(k) >= 0));
            end
        end
    end
    % check that the difference between const(1) - const(2) has the same sign
    % in both tables
    for j = 2:length(t1)
        same(i) = same(i) && (((t1{j-1}.const - t1{j}.const) >= 0) == ((t2{j-1}.const - t2{j}.const) >= 0));
    end
    

    
end



