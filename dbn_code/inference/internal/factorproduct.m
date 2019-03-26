function newf = factorproduct(cpts, varname)
%newf = factorproduct(cpts, varname)
%
% computes the "multiplication" of all factors in CPTS in order to compute
% the JOINT DISTRIBUTION of all conditional probability tables together
%
% INPUT : 
%   CPTS : an array of CondProbTable objects
%   VARNAME : is the string representation of a factor name that is
%       shared by all the input CPTS.
%
% DEBUGGING STATE: should consider replacing by factorproductall() in some
% cases ( where we have many vars to join on per factor pair)
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


%% make sure VARNAME is in each FACTOR
inds = zeros(size(cpts));
for i = 1:length(cpts)
    for j = 1:length(cpts(i).factors)
        if (strcmpi(cpts(i).factors{j}, varname))
            % ...and remember where VARNAME is for later
            inds(i) = j;
            break;
        end
    end
end
temp = cpts(inds > 0);
findices = inds(inds > 0);
    
%% call FACTORPRODUCTWORKER on each pair, repeatedly
if (isempty(temp))
    newf = [];
    return;
end
newfactor = temp(1);
newind = findices(1);
for i = 2:length(findices)
    newfactor = multifactorproductworker(newfactor,temp(i),newind,findices(i));
    % newest factor is always added at the end of the list
    newind = length(newfactor.factors);
end

%% set up return factors
newf = newfactor;
    
    
    
    
    
