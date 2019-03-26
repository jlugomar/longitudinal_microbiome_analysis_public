function b = CheckRunningIntersection(cstree)
% b = CheckRunningIntersection(cstree)
%
% take a cluster set tree and check that it has the running intersection
% property, which is:
% 
% For any two clustersets (i,j), the intersection of i,j is contained in
% each of the clustersets (k) such that i < k < j.
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.

b = zeros(length(cstree));
for i = 1:length(cstree)
    for j = 1:length(cstree)
        if (i == 25 && j == 48)
            foo = true;
        end
        if (i ~= j)
            b(i,j) = dfsIntersectCheck(cstree, i, j);
        end
    end
end
% only zeros should be on the diagonal
zs = find(b == 0);
if (length(zs) > length(cstree))
    error('Running Intersection property is violated after SemiElimTree()');
end

