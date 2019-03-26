function newf = multifactorproductworker(f1,f2, ind1, ind2)
%newf = multifactorproductworker(f1,f2, ind1, ind2)
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

%% find all factors in common: 
% pre-allocate more than we need:
[pairs, inds1, inds2] = intersect(f1.findex, f2.findex, 'legacy');
% make sure the join var is last in these rows:
if (ind1 ~= inds1(end))
    % swap:
    swapinds1 = (ind1 == inds1);
    inds1 = [inds1(~swapinds1), ind1];
    inds2 = [inds2(~swapinds1), ind2];
end

%% find all values of f1.factors(inds1)
allvals = cell(length(inds1),1);
for i = 1:length(inds1)
    allvals{i} = unique([f1.values{inds1(i)}, f2.values{inds2(i)}],'legacy');
end

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

[v, vind, numvals] = ValueIncrement(allvals, []);
for i = 1:numvals    
    if (~skiptable)
        %% this arranges a new conditional probability table
        % need to construct a multi-column match out of all matching pairs:
        shrinktable = f1.table(:,inds1);
        matchtable = shrinktable - repmat(v', size(shrinktable,1),1);
        rows1 = f1.table(sum(matchtable,2) == 0,:);
        p1 = f1.logprob(sum(matchtable,2) == 0);
        shrinktable = f2.table(:,inds2);
        matchtable = shrinktable - repmat(v', size(shrinktable,1),1);
        rows2 = f2.table(sum(matchtable,2) == 0,:);
        p2 = f2.logprob(sum(matchtable,2) == 0);
        r1 = size(rows1,1);
        r2 = size(rows2,1);
        rows1 = repmat(rows1,r2,1);
        rows2 = repmat(rows2,r1,1);
        % complicated indexing achieves staggered repeats for rows1 and rows2
        rrconv = repmat([1:r2:r1*r2]',1,r2) + repmat([0:(r2-1)],r1,1);
        key = reshape(rrconv,size(rows1,1),1);
        rows2 = rows2(key,:);
        index1 = true(1,size(rows1,2));
        index1(inds1) = false;
        index2 = true(1,size(rows2,2));
        index2(inds2) = false;
        % always put the new column last:
        % do this because it is sorted
        newtable = [newtable; rows1(:,index1), rows2(:,index2), rows2(:,inds2)];
    else
        %% if we skipped the table, use computed indices
        index1 = sparsefactorinds(f1, inds1, v);
        index2 = sparsefactorinds(f2, inds2, v);
        p1 = sparse(f1.logprob(index1));
        p2 = sparse(f2.logprob(index2));

        % empty the table:
        newtable = [];
    end
    % multiplication and reshaping actually orders the probabilities in the
    % same way as the index convolution above on the tables.
    %newp = p1' * p2;
    % new logmethod:
    newp = zeros(size(p1));
    srep1 = cell(1,size(p1,2));
    for k = 1:size(p1,2)
        repp1 = repmat(p1(:,k),1,size(p2,2));
        srep1{k} = repp1 + p2;
    end
    for k = 1:size(p1,2)
        for j = 1:size(p2,2)
            newp(k,j) = logplus(srep1{k}(:,j));
        end
    end

    colp = reshape(newp, 1, size(newp,1) * size(newp,2));
    newcolp = [newcolp, colp];
    [v, vind] = ValueIncrement(allvals, v, vind);
end

%% initialize the new factor to have the right values:
newf.table = newtable;
%newf.logprob = newcolp ./ sum(newcolp);
newf.logprob = newcolp;
newf.node = {};

index1 = true(1,length(f1.factors));
index1(inds1) = false;
index2 = true(1,length(f2.factors));
index2(inds2) = false;

newf.factors = {f1.factors{index1}, f2.factors{index2}, f2.factors{inds2}};
newf.values = {f1.values{index1}, f2.values{index2}, f2.values{inds2}};
newf.findex = [f1.findex(index1), f2.findex(index2), f2.findex(inds2)];




    










