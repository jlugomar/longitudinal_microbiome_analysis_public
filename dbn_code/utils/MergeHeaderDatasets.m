function [data,colnames] = MergeHeaderDatasets(d1, cols1, d2, cols2)
%[data,colnames] = MergeHeaderDatasets(d1, cols1, d2, cols2)
%
% Function that merges two datafiles that have some overlapping columns
%
% Returns the NxM data matrix in DATA along with the column labels in a
% cell array of strings COLNAMES.  Only includes columns that appear in
% both datasets.
%
% INPUT:
% D1: dataset one
% COLS1: column labels for dataset one.
% D2: dataset two
% COLS2: column labels for dataset two.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

data = [];
colnames = {};

for i=1:length(cols1)
    c1 = cols1{i};
    c2 = strmatch(c1, cols2, 'exact');
    if (~isempty(c2))
        colnames{end+1} = c1;
        data = [data, [d1(:,i) ; d2(:,c2)]];
    end
end

        
