function newf = factorproductworker(f1,f2, ind1, ind2)
%newf = factorproductworker(f1,f2, ind1, ind2)
%
% computes the joint distribution of two CPTs (F1, F2) that both have the
% same factor at F1.factors{IND1} and F2.factors{IND2}.
%
% INPUT : 
%   F1 : CondProbTable one
%   F2 : CondProbTable two
%   IND1 : indices of F1 that are the same columns as indices of F2 using
%       the index below
%   IND2 : indices of F2
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


MAXTABLE = 100000;

%% error checking:
if (~strcmpi(f1.factors{ind1}, f2.factors{ind2}))
    error('Error, factor mismatch in factorproduct');
end

%% find all values of f1.factors(ind1)
vals = unique([f1.values{ind1}, f2.values{ind2}],'legacy');

newf = f1;

newtable = [];
newcolp = [];


% see if the table exists 
skiptable = false;
if (isempty(f1.table) || isempty(f2.table))
    skiptable = true;
end

% if table is going to be too big, just skip it -
% do indexing via computation instead
if (prod([size(f1.table), size(f2.table)]) > MAXTABLE) 
    skiptable = true;
end

for i = 1:length(vals)    
    v = vals(i);
    if (~skiptable)
        %% this arranges a new conditional probability table
        rows1 = f1.table(f1.table(:,ind1) == v,:);
        rows2 = f2.table(f2.table(:,ind2) == v,:);
        r1 = size(rows1,1);
        r2 = size(rows2,1);
        rows1 = repmat(rows1,r2,1);
        rows2 = repmat(rows2,r1,1);
        % complicated indexing achieves staggered repeats for rows1 and rows2
        rrconv = repmat([1:r2:r1*r2]',1,r2) + repmat([0:(r2-1)],r1,1);
        key = reshape(rrconv,size(rows1,1),1);
        rows2 = rows2(key,:);
        index1 = true(1,size(rows1,2));
        index1(ind1) = false;
        index2 = true(1,size(rows2,2));
        index2(ind2) = false;
        % always put the new column last:
        % do this because it is sorted
        newtable = [newtable; rows1(:,index1), rows2(:,index2), rows2(:,ind2)];

        %% now do the same thing for the actual probabiliteis
        p1 = f1.prob(f1.table(:,ind1) == v);
        p2 = f2.prob(f2.table(:,ind2) == v);
    else
        %% if we skipped the table, use computed indices
        inds1 = sparsefactorinds(f1, ind1, v);
        inds2 = sparsefactorinds(f2, ind2, v);
        p1 = sparse(f1.prob(inds1));
        p2 = sparse(f2.prob(inds2));

        % empty the table:
        newtable = [];
    end
    % multiplication and reshaping actually orders the probabilities in the
    % same way as the index convolution above on the tables.
    newp = p1' * p2;
    colp = reshape(newp, 1, size(newp,1) * size(newp,2));
    newcolp = [newcolp, colp];
end

%% initialize the new factor to have the right values:
newf.table = newtable;
%newf.prob = newcolp ./ sum(newcolp);
newf.prob = newcolp;
newf.node = {};

index1 = true(1,length(f1.factors));
index1(ind1) = false;
index2 = true(1,length(f2.factors));
index2(ind2) = false;

newf.factors = {f1.factors{index1}, f2.factors{index2}, f2.factors{ind2}};
newf.values = {f1.values{index1}, f2.values{index2}, f2.values{ind2}};



    










