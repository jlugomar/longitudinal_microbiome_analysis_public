function [localModel]=learnLocalBN_DiscToDisc(parentData,childData,priorPrecision)
%[localModel]=learnLocalBN_DiscToDisc(parentData,childData,priorPrecision)
%
% This function estimates parameters and computes log-likelihood of the 
% regression-like Bayesian network, given the data. The network assumes a 
% structure where multiple/no discrete parent(s) modulate a single 
% discrete child.
%
% INTERNAL
%
% INPUT:
% PARENTDATA: numeric data matrix of data for the parents of the node
% CHILDDATA: numeric data matrix of data for the one child node
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%
% OUTPUT:
%   localModel.logLLH: the log likelihood of the data, given the child has
%       the input parents
%
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.




%% check input data
numChild = size(childData,1);
numSample = size(childData,2);
numParent = size(parentData,1);
numSampleP = size(parentData,2);

if (numSampleP>0) && (numSample~=numSampleP)
    error('The numbers of samples in parent and child data are not same.');
else
    clear numSampleP;
end

if numChild ~= 1
    error('There must be 1 child node.');
end

%% get configurations of child node
childConfig = unique(childData,'legacy');


%% check priorPrecision
if isempty(priorPrecision) || isempty(priorPrecision.alpha)
	priorPrecision.alpha = 1;
end



%% get configurations of parent nodes
if numParent==0    
    % compute log-likelihood
    % just compute LLH using mmcgeach's method:
    logLLH = LLMultiParent(parentData', childData, priorPrecision.alpha, true);
else
    logLLH = LLMultiParent(parentData', childData, priorPrecision.alpha, false);
end



%% output
localModel.logLLH = logLLH;



return;
