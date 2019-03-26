function [data,colnames] = ReadHeaderDataFile(filename)
% [data,colnames] = ReadHeaderDataFile(filename)
%
% function that reads space-delimited file with a header line of column
% lables.
% returns the NxM data matrix in DATA along with the column labels in a
% cell array of strings COLNAMES
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

%% or, matlab 2013a is terrible at file IO, so use this:
[numdata, cols, strdata, numcolindex, strcolindex] = RCSVLoad(filename, true);
data = numdata;
colnames = cols;

return;


%% read datafile
%data = dlmread(filename,'\t',1,0);

% read new names
%fid = fopen(filename);
% need to skip the %[^\n] construction, which won't work from unix to pc
% and back
%colline = textscan(fid, '%s', 1, 'bufsize', 10000000, 'Delimiter', '\n');
%fclose(fid);
%colnames = textscan(colline{1,1}{1}, '%s', 'bufsize', 10000000);
%colnames = colnames{1};