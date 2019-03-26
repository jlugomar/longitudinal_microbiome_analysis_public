function [localModel]=learnLocalBN_MixToCont(contParentData,discParentData,childData,priorPrecision)
%[localModel]=learnLocalBN_MixToCont(contParentData,discParentData,childData,priorPrecision)
%
% This function estimates parameters and computes log-likelihood of the 
% regression-like Bayesian network, given the data. The network assumes a 
% structure where multiple continuous and discrete parent(s) modulate a 
% single Gaussian child.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.


%% check input arguments
numChild = size(childData,1);
numSample = size(childData,2);
numContParent = size(contParentData,1);
numSampleCP = size(contParentData,2);
numDiscParent = size(discParentData,1);
numSampleDP = size(discParentData,2);

if numChild ~= 1
    error('There must be 1 child node.');
end

if (numSampleDP>0) && (numSample~=numSampleDP)
    error('The numbers of samples in discrete parent and child data are not same.');
else
    clear numSampleDP;
end

if (numSampleCP>0) && (numSample~=numSampleCP)
    error('The numbers of samples in continuous parent and child data are not same.');
else
    clear numSampleCP;
end







%% learn the model parameters
if (numContParent==0) && (numDiscParent==0)
    %% no parent at all
    localModel=learnLocalBN_ContToCont( [], childData, priorPrecision);    
    localModel.discParentConfig=[];
    
elseif (numContParent>0) && (numDiscParent==0)
    %% no discrete parents but continuous parents exist
    localModel=learnLocalBN_ContToCont( contParentData, childData, priorPrecision);    
    localModel.discParentConfig=[];

elseif (numContParent==0) && (numDiscParent>0)
    %% no continuous parents but discrete parents exist
    localModel=learnLocalBN_DiscToCont( discParentData, childData, priorPrecision);    
    
else
    %% both discrete parents and continuous parents exist
    
    %% get discrete parent configurations
    [discParentConfig]=getDiscParentConfig(discParentData);    
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
        tempData = repmat(discParentConfig(:,c),1,numSample);
        tempData = abs(tempData-discParentData);
        tempData = sum(tempData,1);
        tempIdx=find(tempData==0);
        if ~isempty(tempIdx) 
            tempModel=learnLocalBN_ContToCont( contParentData(:,tempIdx), childData(tempIdx), priorPrecision);
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

