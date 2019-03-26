function [discParentConfig]=getDiscParentConfig(parentData)
%[discParentConfig]=getDiscParentConfig(parentData)
%
% 
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.


%% number of parent nodes
numParent = size(parentData,1);


%% values of individual parent nodes
parentValue = cell(1,numParent);


%% cumulative number of parent configurations
numParentConfig = zeros(1,numParent); %preallocate
for p = 1:numParent        
    parentValue{p}=unique(parentData(p,:),'legacy');
    if p==1
        numParentConfig(p) = length(parentValue{p});            
    else
        numParentConfig(p) = numParentConfig(p-1)*length(parentValue{p});
    end
end


%% combine all possible configurations of parent nodes
discParentConfig=nan(numParent,numParentConfig(end)); %preallocate
for p = 1:numParent
    if p==1
        tempRow = parentValue{p};
    else
        tempRow = kron( parentValue{p},ones(1,numParentConfig(p-1)) );            
    end
    discParentConfig(p,:)=repmat(tempRow,1,numParentConfig(end)/numParentConfig(p));
end
clear tempRow;