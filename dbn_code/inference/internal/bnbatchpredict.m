function [pred,correct] = bnbatchpredict_deprc(fname, cpts, varname)
% DEPRICATED
%
% read the datafile FNAME and then predict the value of VARNAME based on
% instantiations of every other variable in the CPTS
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


[data,colnames] = ReadHeaderDataFile(fname);

% find the column corresponding to VARNAME
phencol = 0;
for i=1:length(colnames)
    if (strcmpi(colnames{i}, varname))
        phencol = i;
        break;
    end
end

% find the columns of every variable in the factors in CPTS
[varinds, allvarnames] = GetDataCols(cpts, colnames, varname);


% now treat each line of data as one instantiation
evidence = data(:,varinds);
phn = data(:,phencol);

pred = -1 * ones(size(data,1),1);
for i=1:size(data,1)
    % step 1: reduce all factors, aka, conditional probability tables, to 
    % their evidence
    freds = factorreduce(cpts,allvarnames,evidence(i,:));
   
    % step 2: compute a product of these factors on VARNAME
    fprod = factorproduct(freds, varname);
    
    % step 3: compute the sum of the prob table:
    % predict is the % prob of having the second value of VARNAME
    pred(i) = fprod.prob(2) / sum(fprod.prob);
end

correct = (pred > 0.5) == phn;
        
        
        
        
        
        
        
        
        
        
        
        
        