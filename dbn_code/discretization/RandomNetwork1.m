function [bn, adjmat] = RandomNetwork1(n)
%[bn, adjmat] = RandomNetwork1(n)
%
% generate a specific, arbitrary network, similar to the type we look for
% build them in K2 style: node #n can only have parents numbered 1:(n-1)
% Only netowrk structure is determined here.
%
% useful for illustration purposes.



% use an adjacency matrix
adjmat = zeros(n);
disc = false(1,n);
disc(1:5) = true(1,5);

% create several discrete edges:
adjmat(1,3) = 1;
adjmat(2,3) = 1;
adjmat(1,4) = 1;
adjmat(2,4) = 1;
adjmat(3,4) = 1;
adjmat(4,5) = 1;


% create several gaussian children
c1 = 16;
cm = 20;
adjmat(1,c1:2:cm) = ones(1,(cm - c1)/2 + 1);
adjmat(2,c1:c1+(cm-c1)/2) = ones(1,(cm-c1)/2+1);
adjmat(5,c1:cm) = ones(1,cm - c1 + 1);

% and for those a few more parents!
start = c1;
for c = 6:15
    if (start > cm)
        start = c1;
    end
    adjmat(c,start) = 1;
    start = start + 1;
end

bn = BNfromAdjMat(adjmat,disc);

%% fix up some of the other fields needed in the NODE structure
% for disc nodes, fill in the number of values:
% randomly choose either 2 or three nodes for the number of values
r = rand(1,sum(disc));
r(r > .75) = 3;
r(r < 3) = 2;
% head node should just have 2 values
r(1) = 2;
for i = 1:sum(disc)
    if (r(i) == 2)
        bn(i).values = [1,2];
    end
    if (r(i) == 3)
        bn(i).values = [3,7,11];
    end
end

% do colinds for each node:
for i = 1:n
    bn(i).colind = i;
end

%% write out network to a tgf file:
bayesnet.adjMatrix = adjmat;
nodenames = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};
[MBNodes]=getMarkovBlanketFromBayesNet(bayesnet,nodenames,'1');
outputBayesNetGraph(bayesnet.adjMatrix(MBNodes,MBNodes),nodenames(MBNodes), 'RandomNetwork2');

