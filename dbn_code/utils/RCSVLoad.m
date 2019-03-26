function [numdata, cols, strdata, numcolindex, strcolindex] = ...
    RCSVLoad(fname, usewhitespacedelims, delimiters, FORCEQUOTE, LEADING_NA_ASNUMERIC)
%[numdata, cols, strdata, numcolindex, strcolindex] = RCSVLoad(fname, usewhitespacedelims, delimiters, FORCEQUOTE, LEADING_NA_ASNUMERIC)
%
% slow function for loading csv-files output by the R language.
% These tend to be comma-delimited, but also use quotes to group
% identifiers.  Each column can be a mix of numeric or string data.
%
% INPUT: 
% FNAME: file name to read;
% USEWHITESPACEDELIMS: boolean, will use all whitespace as column delimiters,
%   otherwise will default to using comma as the delimiter
% DELIMITERS: (optional) lets user specify a particular delimiter, like
%   tab.  Should be used with USEWHIESPACEDELIMS = false.
% FORCEQUOTE: (optional) some R csv files will read better if all string
%   columns are treated as quoted string columns.  default = false.
% 
% OUTPUT:
% Output is split into numeric and string data columns:
% NUMDATA : subset of columns in the file that contain only numeric data
% COLS : cell array of strings of numeric column names.
% STRDATA : cell array data for columns containing non-numeric data.
% NUMCOLINDEX : indices of numeric columns.
% STRCOLINDEX : indices of non-numeric columns.
%
% (c) Michael McGeachie 2012.  MIT license. See cgbayesnets_license.txt.

if (nargin <= 1)
    delims = ',';
end

if (nargin >= 2 && usewhitespacedelims)
    delims = [char(32), char(9), char(13)];
else 
    delims = ',';
end

if (nargin >= 3)
    if (strcmp(delimiters, '\t'))
        delims = char(9); %tab
    else
        delims = delimiters;
    end
end

if (nargin < 4)
    FORCEQUOTE = false;
end

if (nargin < 5)
    LEADING_NA_ASNUMERIC = false;
end

fid = fopen(fname);
% need to skip the %[^\n] construction, which is aparently broken in R2013a
header = textscan(fid, '%s', 1, 'Delimiter', '\n');
%header = textscan(fid, '%s', 1, 'bufsize', 10000000, 'Delimiter', '\n');
remain = header{1}{1};

% parse each character in the header
% just use comma and quote as delimiters
done = false;
cols = {};
numempty = 0;
while (~done)
    % use remain(2) here since the standard delimiter is left in the string
    % as the first character.  We skip that one and look at the next
    % character, which might or might not be a doublequote character.
    if (remain(1) == '"' && remain(2) == '"') 
        % empty column name:
        remain = remain(3:end);
        numempty = numempty + 1;
        tok = ['!EMPTY', num2str(numempty),'!'];        
    elseif (remain(2) == '"' && remain(3) == '"')
        % empty column name:
        remain = remain(4:end);
        numempty = numempty + 1;
        tok = ['!EMPTY', num2str(numempty),'!'];
    elseif (remain(1) == '"' || remain(2) == '"')
        [tok, remain] = strtok(remain(2:end),'"');
        % eat the final doublequote character
        remain = remain(2:end);
    else
        [tok, remain] = strtok(remain, delims);
    end
    cols{end+1} = tok;
%    if (isempty(remain) || strcmp('\n',remain(1)))
    if (isempty(remain) || length(remain) < 2)
        done = true;
    end
end

ncols = length(cols);
% parse first line of actual data, after the header
% and make sure the number of columns matches.
% need to skip the %[^\n] construction, which is aparently broken in R2013a
%firstline = textscan(fid, '%s', 1, 'bufsize', 10000000, 'Delimiter', '\n');
firstline = textscan(fid, '%s', 1, 'Delimiter', '\n');
remain = firstline{1}{1};

done = false;
dline = {};
isquoted = false(1,ncols);
i = 1;
while (~done)
    if (remain(1) == '"' && remain(2) == '"') 
        % empty column name:
        isquoted(i) = true;
        remain = remain(3:end);
        tok = '';        
    elseif (remain(2) == '"' && remain(3) == '"')
        % empty column name:
        isquoted(i) = true;
        remain = remain(4:end);
        tok = '';
    elseif (remain(1) == '"' || remain(2) == '"')
        isquoted(i) = true;
        [tok, remain] = strtok(remain(2:end),'"');
        % eat the final doublequote character
        remain = remain(2:end);
    else
        [tok, remain] = strtok(remain, delims);
    end
    dline{end+1} = tok;
    if (isempty(remain) || strcmp('\n',remain(1)))
        done = true;
    end
    i = i + 1;
end
if (length(dline) ~= ncols)
    error('Number of column headers not equal to number of data columns!\n');
end

% guess if each column is going to be numeric or text data
nums = true(1,length(dline));
for i = 1:length(dline)
    numres = str2double(dline{i});
    if (isnan(numres))
        nums(i) = false;
        if (strcmp('NA',dline{i}) && LEADING_NA_ASNUMERIC)
            nums(i) = true;
        end
    end
end
% then build up a parse string for use with textscan()
parsestring = '';
for i = 1:ncols
    if (isquoted(i))
        parsestring = [parsestring, '%q'];
    elseif (nums(i))
        parsestring = [parsestring, '%f'];
    else
        if (FORCEQUOTE)
            parsestring = [parsestring, '%q'];
        else
            parsestring = [parsestring, '%s'];
        end
    end        
end
parsestring = [parsestring];
%[tdata, pos] = textscan(fid, parsestring, 'Delimiter', delims, 'TreatAsEmpty', {'NA','""'}, ...
%    'BufSize', 50000000, 'MultipleDelimsAsOne', 0);
[tdata, pos] = textscan(fid, parsestring, 'Delimiter', delims, 'TreatAsEmpty', {'NA','""'}, ...
    'MultipleDelimsAsOne', 0);

% check we got the same number of rows in each column
prevs = length(tdata{1});
for i = 1:length(tdata)
    s = length(tdata{i});
    if (s ~= prevs)
        error('Malformed R file at position %d in data column %d, row %d\n', pos, i, s+1);
    end
    prevs = s;
end
numdata = zeros(length(tdata{1})+1,sum(nums));
strdata = cell(1,ncols - sum(nums));
dcol = 1;
scol = 1;
numcolindex = [];
strcolindex = [];
% assemble output
for i = 1:ncols
    if (nums(i))
        numdata(1,dcol) = str2double(dline{i});
        try
            numdata(2:end,dcol) = tdata{i};
            numcolindex(dcol) = i;
            dcol = dcol + 1;
        catch err
            errmsg = err.message;
            if (strfind(errmsg, 'Conversion to double'))
                % error converting this column to numeric data; treat it as
                % string data instead:
                strdata{scol} = {dline{i}, tdata{i}{:}}';
                strcolindex(scol) = i;
                scol = scol + 1;
                % then remove this from the numeric data:
                numdata = numdata(:,[1:dcol-1,dcol+1:end]);
            end
        end
    else
        strdata{scol} = {dline{i},tdata{i}{:}}';
        strcolindex(scol) = i;
        scol = scol + 1;
    end
end


fclose(fid);


