function [data, cols] = LoadCSVFile(fname, delim)
%[data, cols] = LoadCSVFile(fname, delim)
%
% function to load a CSV or similar file, and grab both the column headers
% and the numeric data columsn. 
%
% INPUT: 
% FNAME: file name to read from 
% DELIM: delimter, default is the comma (',') other suggestions might be
%   the tab '\t'
%
% OUTPUT:
% DATA: data columns from the file
% COLS: column headers; variables names.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    [numdata, cols, strdata, numcolindex, strcolindex] = RCSVLoad(fname, false, ',');
end
if (nargin < 3)
    [numdata, cols, strdata, numcolindex, strcolindex] = RCSVLoad(fname, false, delim);
end

data = numdata;
return;

