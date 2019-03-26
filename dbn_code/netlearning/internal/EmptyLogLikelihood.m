function LLH = EmptyLogLikelihood(data, disc, priorPrecision)
% function LLH = EmptyLogLikelihood(data, disc)
%
% computes the base loglikelihood of an empty network over a dataset
% 
% INPUT:
%   DATA: the dataset itself
%   DISC: boolean indicator for each column of DATA; DISC(i) is true if
%       DATA(:,i) is a discrete variable.
%
% OUTPUT:
%   LLH: loglikelihood of the data
%
%
%
% Copyright Michael McGeachie, 2015.  MIT license. See cgbayesnets_license.txt.



% loop through each variable

LLH = 0;
ncols = size(data,2);
for i = 1:ncols    
    % call bayesfactor_worker()
    [bf, basellh] = bayesfactor_worker([],[],data(:,i),[],[],disc(i),priorPrecision);
    LLH = LLH + basellh;
end






