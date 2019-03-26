function nodes = BNfromAdjMat(mat, discrete, nodenames)
%nodes = BNfromAdjMat(mat, discrete, nodenames)
%
% Takes as input an adjacency matrix and output the array of NODE structures
% that corresponds to a hybrid discrete/continuous Bayesian Network with
% the same connectivity.
% Maintains order of nodes according to listing in MAT.
%
% INPUT : 
%   MAT : adjacency matrix, such that if MAT(i,j) is non-zero, there is a
%       link from NODES(i) to NODES(j)
%   DISCRETE : boolean array indicating which nodes are discrete
%   NODENAMES : optional names for each of the nodes
%
% OUTPUT : 
%   NODES : array of NODE structures that encodes the Bayes Net.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    discrete = true(1, size(mat,1));
end

if (nargin < 3)
    nodenames = cell(size(discrete));
    for i=1:length(discrete)
        nodenames{i} = num2str(i);
    end
end

nodes = repmat(Node(),1,size(mat,1));

for i = 1:length(nodes)
    nodes(i).self = nodenames{i};
    nodes(i).discrete = discrete(i);
end


% should also work with sparse matrices:
[rowinds,colinds] = find(mat);
for i = 1:length(rowinds)
    if (isempty(nodes(rowinds(i)).children))
        nodes(rowinds(i)).children = {nodes(colinds(i)).self};
        nodes(rowinds(i)).cindex = colinds(i);
    else
        nodes(rowinds(i)).children{end+1} = nodes(colinds(i)).self;
        nodes(rowinds(i)).cindex(end+1) = colinds(i);
    end

    if (isempty(nodes(colinds(i)).parents))
        nodes(colinds(i)).parents = {nodes(rowinds(i)).self};
        nodes(colinds(i)).pindex = rowinds(i);
     else
        nodes(colinds(i)).parents = {nodes(colinds(i)).parents{:},nodes(rowinds(i)).self};
        nodes(colinds(i)).pindex(end+1) = rowinds(i);
    end
end

