function [varinds, allvarnames] = GetDataCols(cpts, colnames, varname)
%[varinds, allvarnames] = GetDataCols(cpts, colnames, varname)
%
% find the columns of every variable in the factors in CPTS
%
% INPUT : 
%   CPTS : array of CondProbTables that represent each node in the Bayesian
%       network
%   COLNAMES : cell array of strings that are the column headers for data
%       and the names of each variable
%   VARNAME : the name of the phenotype var that is being predicted
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


allvarnames = {};
nextind = 1;
for i=1:length(cpts)
    c = cpts(i);
    for j=1:length(c.factors)
        name = c.factors{j};
        % make sure this factor is NOT varname
        if (strcmpi(name, varname))
            continue;
        end
        % check if we've seen this factor before
        found = false;
        for k=1:length(allvarnames)
            if (strcmpi(allvarnames{k},name))
                found = true;
                break;
            end
        end
        if (~found)
            allvarnames{nextind} = name;
            nextind = nextind + 1;
        end
    end
end
varinds = zeros(size(allvarnames));
for i=1:length(allvarnames)
    for j=1:length(colnames)
        if (strcmpi(colnames{j}, allvarnames{i}))
            varinds(i) = j;
            break;
        end
    end
end
