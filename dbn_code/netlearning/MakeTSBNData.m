function [d0,dn,nsids] = MakeTSBNData(data, subjectids)
% tsdata = MakeTSBNData(data, subjectids)
% 
% assemble a 2-stage Dynamic BN dataset from times series data
%
% takes input data and arranges it so the first time a subject id is
% encountered, its slotted into the T0 data.  The second time its
% encountered, its slotted into both the TN and the T0 data.  The last time
% its encountered, its only slotted into the TN data.
%
% INPUT :
%   DATA : N x M datamatrix with M vars in columns and N observations in rows
%   SUBJECTIDS : N x 1 array of subject ids, for matching rows with 
%       observationss from the same subject at different time points. 
%       Parallel to DATA.
%   TIMES : N x 1 array of integeres indicating observations at different
%       time slices. Parallel to DATA.
%
% OUTPUT :
%   d0 : all data except the last time point for each subject.
%   dn : all data except the first time point for each subject.
%   nsids : parallel array of subject ids for the output data.  The k'th
%       subject's data are in d0(nsids == k) and dn(nsids == k).
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

% this also sorts from lowest to highest
subs = unique(subjectids);

% take input data and arrange it so the first time a subject id is
% encountered, its slotted into the T0 data.  The second time its
% encountered, its slotted into both the TN and the T0 data.  The last time
% its encountered, its only slotted into the TN data.
d0 = [];
dn = [];
nsids = [];
for i = 1:length(subs)
    si = subjectids == subs(i);
    ds = data(si,:);
    % drop any subject that only appears in one time point
    if (sum(si) > 1)
        sarray = subs(i) * ones(size(ds(:,1)));
        d0 = [d0;ds(1:end-1,:)];
        dn = [dn;ds(2:end,:)];
        nsids = [nsids;sarray(2:end)];
    end
end




