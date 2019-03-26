function MB = bdnmarkovblanket(node, nodename)
% MB = bdnmarkovblanket(node, nodename)
%   find the markov blanket (MB) of node named NODENAME
%   MB = parents + children + other parents of children
%   NOTE: case-insensitive 


namednode = bdnfind(node,nodename);
if isempty(namednode)
    error(['node "' nodename '" not found in: ' node(1).filename]);
end

% add the children of the node
MB = namednode.children;

% add the other parents of children of the node
children = namednode.children;
nchildren = size(children,2);
for i=1:nchildren;
    cnode = bdnfind(node, children{i});
    % debug cnode.self
    ncparents = length(cnode.parents);
    for j=1:ncparents
        if strcmpi(cnode.parents{j}, 'NIL')
            % there is no parent, do nothing
        elseif strcmpi(cnode.parents{j},nodename)
            % its the named node, so do nothing
        else
            % its an other parent, add to the list
            if isempty(strmatch(cnode.parents{j}',MB))
                MB{size(MB,2)+1} = cnode.parents{j};
            end
        end
        % debug MB
    end
end

% add the parents of the node

nparents = length(namednode.parents);
for i=1:nparents
    if strcmpi(namednode.parents{i}, 'NIL')
        % its the root, do nothing
    else
        % TK add namednode.parents to the MB };
        if isempty(strmatch(namednode.parents{i}',MB))
            MB{size(MB,2)+1} = namednode.parents{i};
        end
    end
end

function namednode = bdnfind(node, nodename)

% find node named nodename

namednode = []; % default value if not found

nodes=length(node);

for i=1:nodes
    % find node named nodename
    if strcmpi(node(i).self, nodename) % allow upper/lower mismatch
        namednode=node(i);
        return
    end
end