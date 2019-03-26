function [bf,independence_bf] = bayesfactor_worker(parentdata, newparentdata, childdata, ...
    discparent, discnewparent, discchild, priorPrecision)
%bf = bayesfactor_worker(parentdata, newparentdata, childdata, ...
%    discparent, discnewparent, discchild, priorPrecision)
% 
% take various parent configurations for a data column and figure out which
% likelihood function to use to compute the likelihood.
%
% INPUT
%   PARENTDATA: a matrix of the data for the existing parents of the child.
%       should be in a [samples x vars] matrix.  Can be empty.
%   NEWPARENTDATA: a matrix of the data for the new parents of the child.
%       should be in a [samples x vars] matrix.  Can be empty.
%   CHILDDATA: a column of data for the child. 
%   DISCPARENT: a boolean array parallel to PARENTDATA, indicating if each
%       column in PARENT data is discrete or not.
%   DISCNEWPARENT: a boolean array parallel to NEWPARENTDATA, indicating if each
%       column in NEWPARENT data is discrete or not.
%   DISCCHILD: a boolean value indicating if CHILDDATA is discrete or not.
%   PRIORPRECISION: a structure including the usual HybridBayesNets
%     parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%
% OUTPUT: 
%   BF: number indicating the difference in log likelihood for the
%       connections indicated.
%
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

if (~isempty(parentdata))
    cparent = parentdata(:,~discparent)';
    dparent = parentdata(:,discparent)';
else
    cparent = [];
    dparent = [];
end
if (~isempty(newparentdata))
    ncparent = [cparent; newparentdata(:,~discnewparent)'];
    ndparent = [dparent; newparentdata(:,discnewparent)'];
else
    ncparent = cparent;
    ndparent = dparent;
end


if (~discchild)
    % use mixed to continuous likelihood computation
    [tempModel]=learnLocalBN_MixToCont(cparent,dparent,childdata',priorPrecision);
    tempLogLLH_Ind = tempModel.logLLH;
    [tempModel]=learnLocalBN_MixToCont(ncparent,ndparent,childdata',priorPrecision);
    tempLogLLH_Dep = tempModel.logLLH;
else
    % use discrete to discrete likelihood computation
    [tempModel]=learnLocalBN_DiscToDisc(dparent, childdata', priorPrecision);    
    tempLogLLH_Ind = tempModel.logLLH;
    [tempModel]=learnLocalBN_DiscToDisc(ndparent, childdata', priorPrecision);    
    tempLogLLH_Dep = tempModel.logLLH;    
end
bf=tempLogLLH_Dep-tempLogLLH_Ind;
independence_bf=tempLogLLH_Ind;

