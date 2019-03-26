function OutputForWEKA(fname, data, cols, uniquevals, disc)
%OutputForWeka(fname, data, cols, uniquevals, disc)
%
% Writes a dataset to an ARFF file that is readable by Weka 3.6.9. 
%
%
% INPUT
% FNAME: output file name.  Should end in ".arff"
% DATA: data matrix, variables stored in columns
% COLS: column names for the variables in DATA.
% UNIQUEVALS: discrete nodes must come with a list of unique values they
%   obtain.
% DISC: boolean array parallel to NAMES indicating if the data column is
%   discrete or continuous
%
% (c) Michael McGeachie, 2013, all rights reserved.

fid = fopen(fname, 'w');

% name the file/data/experiment
fprintf(fid, '@relation %s\n', fname);

% print out info on each variable
aorder = [2:length(cols),1];
for i = aorder
    fprintf(fid, '@attribute %s ', cols{i});
    if (disc(i))
        uvals = uniquevals{i};
        fprintf(fid, '{');
        for j = 1:length(uvals)
            % may need to tweak this so training and test datasets come out
            % with the same list of unique values
            fprintf(fid, '%s ', num2str(uvals(j)));
        end
        fprintf(fid, '}\n');
    else
        fprintf(fid, 'numeric\n');
    end
end

% now print out the main data, separated by commas, but no spaces
fprintf(fid,'@data\n');
for j = 1:size(data,1)
    for i = 1:length(aorder)
        fprintf(fid,'%s', num2str(data(j,aorder(i))));
        if i < length(aorder)
            fprintf(fid,',');
        end
    end
    fprintf(fid,'\n');
end
