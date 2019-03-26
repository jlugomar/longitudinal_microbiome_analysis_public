function [BN] = LearnParamsBN(BN)
%[BN] = LearnParams(BN)
%
% Main function for learning the parameters of a CG Bayesian Network. Use
% after learning the structure of the Bayesian Network.
% 
%  INPUT: 
%  BN: The Bayesian Network to learn parameters for.  Class "BayesNet"
%   object.
%
%  OUTPUT:
%  BN: the BayesNet class object is updated to have meaningful fields in:
%     BN.tree 
%     BN.nodes
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.



[treeOut, nodes] = LearnParams(BN.nodes, '', BN.data, BN.cols, BN.priorPrecision, ...
    BN.discvals, BN.disc);

BN.tree = treeOut;
BN.nodes = nodes;
