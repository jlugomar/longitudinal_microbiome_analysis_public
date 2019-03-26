function newf = joindisjointcpts(f1,f2)
%newf = joindisjointcpts(f1,f2)
%
% computes the joint distribution of two CPTs (F1, F2) that have an empty 
% intersection.
% INPUT : 
%   F1 : CondProbTable class object
%   F2 : CondProbTable class object
%
% also just returns the superset CTP in the case of one the inputs being a
% subset of the other.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


MAXTABLE = 100000;
newf = f2;

newtable = [];
newcolp = [];

% see if these factors are actually disjoint:
if (isempty(setdiff(f1.findex, f2.findex, 'legacy')))
    % everything in f1 is contained in f2:
    newf = f2;
    return;
elseif (isempty(setdiff(f2.findex, f1.findex, 'legacy')))
    newf = f1;
    return;
end

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

if (~skiptable)
    %% this arranges a new conditional probability table
    rows1 = f1.table;
    rows2 = f2.table;
    r1 = size(rows1,1);
    r2 = size(rows2,1);
    rows1 = repmat(rows1,r2,1);
    rows2 = repmat(rows2,r1,1);
    % complicated indexing achieves staggered repeats for rows1 and rows2
    rrconv = repmat([1:r2:r1*r2]',1,r2) + repmat([0:(r2-1)],r1,1);
    key = reshape(rrconv,size(rows1,1),1);
    % just re-order the rows of table 2
    rows2 = rows2(key,:);
    % always put the new column last:
    % do this because it is sorted
    newtable = [rows1, rows2];
%else
%   error(['Conditional Prob Table Overflow!  Size > ', num2str(MAXTABLE)]);
end
% select the probabilities corresponding to the rows 
p1 = f1.logprob;
p2 = f2.logprob;
% multiplication and reshaping actually orders the probabilities in the
% same way as the index convolution above on the tables.
%newp = p1' * p2;
% new logmethod:
newp = zeros(size(p1));
srep1 = cell(1,size(p1,2));
for i = 1:size(p1,2)
    repp1 = repmat(p1(:,i),1,size(p2,2));
    srep1{i} = repp1 + p2;
end
for i = 1:size(p1,2)
    for j = 1:size(p2,2)
        newp(i,j) = logplus(srep1{i}(:,j));
    end
end


colp = reshape(newp, 1, size(newp,1) * size(newp,2));
newcolp = [newcolp, colp];

%% initialize the new factor to have the right values:
newf.table = newtable;
%newf.logprob = newcolp ./ sum(newcolp);
newf.logprob = newcolp;
newf.node = [];

newf.factors = {f1.factors{:}, f2.factors{:}};
newf.values = {f1.values{:}, f2.values{:}};
newf.findex = [f1.findex, f2.findex];



