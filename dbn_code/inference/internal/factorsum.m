function newf = factorsum(cpt, varname)
%
% computes the "summation" of a factor in CPT 
%
% INPUT : 
%   CPT : a ConditProbTable
%   VARNAME : string representation of a factor name in CPT.
%
% When a factor has a variable "summed out" it is removed from the
% conditional probability table by combining all elements of the table that
% have the same assignments to the other vars (in CPT.FACTORS - VARNAME).
%
% Summing and Multiplying factors are the two basic operations of the
% variable elimination algorithm for inference in bayesian networks.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


%% make sure VARNAME is in the FACTOR
varind = 0;
for i = 1:length(cpt.factors)
    if (strcmpi(cpt.factors{i}, varname))
        varind = i;
        break;
    end
end

%% sum the cpt.prob elements for matching values of (factors - VARNAME)
vals = cpt.values{varind};
finds = true(1,length(cpt.factors));
finds(varind) = false;

% easy to do using our SKIP-TABLES computational indexing
for i = 1:length(vals)
    v = vals(i);
    inds = sparsefactorinds(cpt, varind, v);
    if (i > 1)
        newprob = newprob + sparse(cpt.prob(inds));
    else
        newprob = sparse(cpt.prob(inds));
    end
end

%% build the new table for the new factor
newtable = [];
if (~isempty(cpt.table))
    % in the KEEP-TABLES case, just keep the rows corresponding to a single
    % value of VARNAME, delete the rest.  Also delete the column VARNAME
    v = vals(1);
    newtable = cpt.table(cpt.table(:,varind) == v,:);
    newtable = newtable(:,finds);
end

%% set up new CPT:
newf = CondProbTable(cpt.node, {cpt.factors{finds}}, newtable, newprob, {cpt.values{finds}});




















