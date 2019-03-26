function [bn, adjmat] = RandNet(n, ndisc)
%[bn, adjmat] = RandNet(n)
%
% Generates a random bayesian network, with some discrete and continuous
% nodes. Networks are built so that node #n can only have parents numbered 
% 1:(n-1), meaning the nodes of the BN are already in topological order.
%
% INPUT
% N: number of nodes in the network
% NDISC: optional: number of discrete nodes, out of the total N.  The
%   N-NDISC remaining nodes will be Gaussian.
%
% OUTPUT
% BN: array of NODE objects representing the Bayes Net
% ADJMAT: adjacency matrix of the network generated.
%
% (c) Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.


DISC_EDGE_RATE = .4;
CONT_EDGE_RATE = .25;

% use an adjacency matrix
adjmat = zeros(n);
disc = false(1,n);



if (nargin < 2)
    ndisc = randn(1,1) * n/5 + n *.4;
    ndisc = floor(ndisc);
    ndisc = max(ndisc,1);
    ndisc = min(ndisc,n-5);
end

% create first discrete nodes
disc(1:ndisc) = true;

% create a few discrete node edges
for i = 1:ndisc
    r = rand(1,ndisc);
    r(1:i) = ones(1,i);
    adjmat(i, r <= DISC_EDGE_RATE) = 1;
end

% deal with gaussian nodes
for i = 1:n
    r = rand(1,n);
    r(1:i) = ones(1,i);
    adjmat(i, r <= CONT_EDGE_RATE) = 1;
end

bn = BNfromAdjMat(adjmat,disc);

%% fix up some of the other fields needed in the NODE structure
% for disc nodes, fill in the number of values:
% randomly choose either 2 or three nodes for the number of values
r = rand(1,ndisc);
r(r > .66) = 3;
r(r < 3) = 2;
% head node should just have 2 values
r(1) = 2;
for i = 1:ndisc
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

