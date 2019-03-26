function parent = addedge(parent,child)
% parent = addedge(parent,child)
% 
% helper function for various network operations upon NODES
%
% INPUT : 
%   PARENT : node struct
%   CHILD : node struct
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

flag = false;
if (strcmpi(parent.self, child.self))
    flag = true;
    return;
end
for i = 1:length(parent.children)
    if (strcmpi(parent.children{i}, child.self))
        flag = true;
        break;
    end
end
if (~flag)
    parent.children{end+1} = child.self;
    if(isempty(parent.cindex))
        parent.cindex = child.index;
    else
        parent.cindex(end+1) = child.index;
    end
end