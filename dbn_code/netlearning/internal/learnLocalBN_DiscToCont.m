function [localModel]=learnLocalBN_DiscToCont(parentData,childData,priorPrecision)
%[localModel]=learnLocalBN_DiscToCont(parentData,childData,priorPrecision)
%
% This function estimates parameters and computes log-likelihood of the 
% regression-like Bayesian network, given the data. The network assumes a 
% structure where multiple/no discrete parent(s) modulate a single 
% Gaussian child.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.


%% check input arguments
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



%% learn the model parameters
if numParent ==0    
    %% no parents at all
    localModel=learnLocalBN_ContToCont( [], childData, priorPrecision);    
    localModel.discParentConfig=[];

else

    %% get discrete parent configurations
    [discParentConfig]=getDiscParentConfig(parentData);    
    numParentConfig = size(discParentConfig,2);
    
    %% update hyperparemeter "nu" in priorPrecision
    priorPrecision.nu = priorPrecision.nu/numParentConfig;
    
    %% go over all the configurations of the discrete parent. 
    postMean = cell(1,numParentConfig);
    postVar = cell(1,numParentConfig);
    postPrec = cell(1,numParentConfig);    
    regreCoeff = cell(1,numParentConfig);
    logLLH = zeros(1,numParentConfig);
    for c = 1:numParentConfig
        tempData = repmat(discParentConfig(:,c),1,numSample );
        tempData = abs(tempData-parentData);
        tempData = sum(tempData,1);
        tempIdx=find(tempData==0);
        if ~isempty(tempIdx) 
            tempModel=learnLocalBN_ContToCont( [], childData(tempIdx), priorPrecision);
            postMean{c} = tempModel.postMean;
            postVar{c} = tempModel.postVar;
            postPrec{c} = tempModel.postPrec;
            regreCoeff{c} = tempModel.regreCoeff;    
            logLLH(c) = tempModel.logLLH;
        else
            postMean{c} = 1;
            postVar{c} = 1;
            postPrec{c} = 1;
            regreCoeff{c} = 1;    
            logLLH(c) = 0;
        end
        clear temp*;
    end
    
    
    %% output
    localModel.discParentConfig=discParentConfig;
    localModel.postMean = postMean;
    localModel.postVar = postVar;
    localModel.postPrec = postPrec; 
    localModel.regreCoeff = regreCoeff;
    localModel.logLLH = sum(logLLH);
    
end





return;

