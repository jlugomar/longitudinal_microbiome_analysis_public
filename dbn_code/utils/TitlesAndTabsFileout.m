function TitlesAndTabsFileout(fname, titlearray, data, precision)
%TitlesAndTabsFileout(fname, titlearray, data)
%
% function to write out a data array along with column headers for those
% data columns.
%
% INPUT: 
% FNAME: file name to write to
% TITLEARRAY: column headers, cell array of strings
% DATA: data array, 1 variable per column.  Numeric data only.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 4)
    precision = 10;
end
if (isempty(data))
    return;
end

fid = fopen(fname,'wt');
for i = 1:length(titlearray)
    fprintf(fid, titlearray{i});
    if (i ~= length(titlearray))
        fprintf(fid, '\t');
    end 
end
fprintf(fid, '\n');
fclose(fid);
dlmwrite(fname,data,'-append','delimiter', '\t', 'newline','pc','precision', precision);