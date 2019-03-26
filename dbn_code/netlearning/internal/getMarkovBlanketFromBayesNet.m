function [MBNodes]=getMarkovBlanketFromBayesNet(BayesNet,nodeNames,targetNodeName)
%[MBNodes]=getMarkovBlanketFromBayesNet(BayesNet,nodeNames,targetNodeName)
%
% Computes a Markov Blanket for a given CG Bayesian network.  This is all
% the nodes that are required to predict the value of TARGETNODENAME.
%
% INPUT:
% BAYESNET: structure with:
%   BayesNet.adjMatrix: an adjacency matrix representing the Bayes Net.
% NODENAMES: cell array of strings listing the node names
% TARGETNODENAME: name of the node which should be the center of the Markov
%   Blanket.
%
% OUTPUT:
% MBNODES: List of node indices comprising the Markov Blanket.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.


%% index of target node
targetNode = strmatch(targetNodeName, nodeNames, 'exact');


%% get nodes under Markov blanket
parentNodes = find(BayesNet.adjMatrix(:,targetNode));
childNodes = find(BayesNet.adjMatrix(targetNode,:));
parentOfChildNodes = [];
for c = childNodes
    parentOfChildNodes = union(parentOfChildNodes, find(BayesNet.adjMatrix(:,c)),'legacy');    
end


%% get all nodes under the Markov blanket
MBNodes = union(parentNodes,childNodes,'legacy');
MBNodes = union(MBNodes,parentOfChildNodes,'legacy');

