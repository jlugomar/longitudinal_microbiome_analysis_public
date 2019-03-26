function cpt = InitDiscreteFactor(cst,node,nodes,data, priorPrecision)
%cpt = InitDiscreteFactor(cst,node,nodes,data, priorPrecision)
%
% this function works for initializing discrete nodes in either HYBRID or
% DISCRETE bayesian networks.  for discrete networks, call this with 
% CST = [].
%
% INPUT : 
%   CST : a ClusterSetTree class object.
%   NODE : one of the NODES
%   NODES : array of class NODE instances, describing the Bayesian Netowrk
%   DATA : array of data to initialize the CondProbTable
%   PRIORPRECISION : priorPrecision.alpha is the psuedo counts of prior
%       data seen.
%
%
% OUTPUT : CPT is an instantiation of class CondProbTable, representing the
% conditional probability table for a discrete node in a bayesian network
% or a discrete ClusterSetTree in a hybrid bayesian network's junction tree
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


cpt = CondProbTable();

cpt.node = node.self;
cpt.factors = {node.self};
cpt.findex = [node.index];
if (~isempty(cst))
    for i = 1:length(cst.members)
        % 5/11/10: changed this to remove the {} on the added factor:
        % also check to make sure we haven't already added this node,
        % if it is identical with the input NODE
        if (~strcmpi(node.self, nodes(cst.members(i)).self))
            cpt.factors{end+1} = nodes(cst.members(i)).self;
            cpt.findex(end+1) = cst.members(i);
        end
    end
else
    cpt.factors = {node.self, node.parents{1:end}};
end

%% find indices to factor members in the datafile/dataset:
cols = zeros(1,length(nodes));
for i = 1:length(cpt.factors)
    cols(i) = nodes(cpt.findex(i)).colind;
end
% remove unfound columns:
cols = cols(cols > 0);
d = data(:,cols);

%% build table:
t = [];
for i = 1:length(cpt.findex)
    cpt.values{i} = cell2mat(nodes(cpt.findex(i)).values);
    t = repmat(t,length(cpt.values{i}),1);
    % make sure we get at least one copy of the values :
    r = max(1,size(t,1) / length(cpt.values{i})); 
    newv = repmat(cpt.values{i}, 1, r);
    newcol = reshape(newv',length(cpt.values{i}) * r,1);
    t = [t, newcol];
end
cpt.table = t;
if (isempty(cpt.logprob))
    cpt.logprob = zeros(size(cpt.table,1),1);
end

%% sum up number of co-occurrences of TABLE and put them in LOGPROB
% init prior weights on every table entry: 
cpt.logprob = cpt.logprob + priorPrecision.alpha / (size(cpt.table,1) * size(cpt.table,2));
if (size(d,1) > size(cpt.table,1))    
    for i = 1:size(cpt.table,1)
        assignment = cpt.table(i,:);
        a = repmat(assignment,size(d,1),1);
        matches = d - a;
        % matches are those rows identically zero
        zs = sum(abs(matches),2);
        % so we can sum the zeros and count the matches
        cpt.logprob(i) = cpt.logprob(i) + sum(zs == 0);
    end
else
    for i = 1:size(d,1)
        % just increment each prob for every datapoint
        assignment = d(i,:);
        a = repmat(assignment,size(cpt.table,1),1);
        matches = cpt.table - a;
        % the match is the row that's zero
        zs = logical(sum(abs(matches),2) == 0);
        % increment that one
        cpt.logprob(zs) = cpt.logprob(zs) + 1;
    end    
end

%% generally don't bother with normalizing anything
%a = sum(cpt.logprob);
%cpt.logprob = cpt.logprob ./ a;
v = length(cpt.values{1});
for i = 1:size(cpt.table,1)/v
    % get every v^th + i row, for v = length(cpt.values{1})
    inds = [((i-1)*v)+1:i*v];
    a = sum(cpt.logprob(inds));
    if (a ~= 0)
        cpt.logprob(inds) = cpt.logprob(inds) ./ a;
    end
end
% finally take log of all probs:
if (cpt.logprob <= 0)
    error('Data too fragmented!  Must increase prior to avoid zero probability occurrences!');
end
% on the other hand, also an error if cpt.logprob is > 1 :
if (cpt.logprob > 1)
    error('Raw prob > 1, expecting probability in range [0,1]');
end
cpt.logprob = log(cpt.logprob);




