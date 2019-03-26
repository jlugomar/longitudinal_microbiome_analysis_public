function nodes = nodetree(nodes, data, colnames, discvals, is_discrete)
%nodes = nodetree(nodes, data, colnames, discvals, is_discrete)
%
% takes an array of NODES and add indexing arrays to their parents,
% children, and self.
%
% also remove the 'NIL' parent notifier for roots
%
% If inputs 2 and 3 are included, this will attempt to determine if each
% variable is DISCRETE or not, and if discrete, what are the VALUES of that
% variable
% 
% input : 
%   NODES : array of class NODE instantes 
%   DATA : data for inference of teh bayes net parameters
%   COLNAMES : names of cols, matching NODES().self
%   DISCVALS : (optional) cell array of values for each discrete node,
%       incase it is helpful to specify all possible values for nodes that
%       have extremely rare values, which might not exist in all datasets.
%       Parallel to COLNAMES
%   IS_DISCRETE : (optional) logical array parallel to COLNAMES indicating
%       which variables are discrete
%
% FILLS IN THESE FIELDS OF NODES():
% node.index = integer index into nodes array for itself
% node.pindex = array of indices into nodes array for parents
% node.cindex = array of indices into nodes array for children
% node.discrete = TRUE if DATA is input and can be discerned from the DATA
% node.colind = index into COLNAMES such that node(i).self = colnames(node(i).colind)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

% skip everything related to data if those arguments are not included:
if (nargin < 2)
    skipdata = true;
else
    skipdata = false;
end

% if we have column names, make sure that they correspond to the node names
if (nargin >= 3)
    nodes = AddColIndstoNodes(nodes, colnames);
    if (isempty(data))
        skipdata = true;
    end
end

% first set the self index:
for i=1:length(nodes)
    nodes(i).index = i;
end

if (nargin < 5)
    is_discrete = [];
end

% convert from string pointers to other nodes to node index pointers
for i=1:length(nodes)
    if (~isempty(nodes(i).parents))
        if (strcmpi(nodes(i).parents{1}, 'NIL'))
            nodes(i).parents = {};
            nodes(i).pindex = [];
        else
            for p=1:length(nodes(i).parents)
                % if the parent index list isn't already populated, add the
                % parents:
                if (isempty(nodes(i).pindex))
                    nodes(i).pindex = nodeindex(nodes,nodes(i).parents{p});
                elseif (isempty(find(nodes(i).pindex == nodeindex(nodes,nodes(i).parents{p}),1)))
                    nodes(i).pindex(end+1) = nodeindex(nodes,nodes(i).parents{p});
                end
            end
        end
    end
    nodes(i).cindex = [];
    for c=1:length(nodes(i).children)
        nodes(i).cindex(end+1) = nodeindex(nodes,nodes(i).children{c});
    end
end

%% find column name indices
for i = 1:length(nodes)
    found = false;
    for j = 1:length(colnames)
        if (strcmpi(colnames{j}, nodes(i).self))
            nodes(i).colind = j;
            found = true;
            break;
        end
    end
    if (~found)
        nodes(i).colind = -1;
    end
end

%% give default values if we don't have data from which to compute these:
if (skipdata)
    %nodes(i).discrete = [];
    nodes(i).colind = -1;
else %if (~skipdata)
    %% figure out if each variable is discrete:
    isdisc = IsDiscrete(data);
    for i = 1:length(nodes)
        if (isempty(is_discrete))
            nodes(i).discrete = isdisc(nodes(i).colind);
        else
            nodes(i).discrete = is_discrete(nodes(i).colind);
        end
    end

    %% count values for each discrete variable:
    for i = 1:length(nodes)
        if (nodes(i).discrete)
            nodes(i).values = num2cell(unique(data(:, nodes(i).colind),'legacy'));
            if (nargin >=4)
                if (~isempty(discvals) && ~isempty(discvals{nodes(i).colind}))
                    if (length(discvals{nodes(i).colind}) > length(nodes(i).values))
                        nodes(i).values = discvals{nodes(i).colind};
                    end
                end
            end
        end
    end

end



end

%% node indexing:
% helper function that gets the index of a node into the tree by searching
% on the name of the node:
function index = nodeindex(nodes, nodename)
for index = 1:length(nodes)
    if (strcmpi(nodes(index).self, nodename))
        return;
    end
end
index = 0;       
end