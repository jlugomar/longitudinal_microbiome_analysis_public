function [rint, hascycle] = dfsIntersectCheck(tree, start, stop)
% Perform a depth-first search of the graph starting from 'start'.
% [rint, hascycle] = dfsIntersectCheck(adj_mat, start, stop)
%
% Input:
% 'tree' is a CLUSTERSETTREE array
% 'start' is the index into tree to start the search 
% 'stop' is the index into tree to stop the search
%
% Output:
% 'hascycle' is true iff a directed cycle is found.
% 'rint' is 1 if the running intersection property holds from tree nodes
%       TREE(START) to TREE(STOP).  -1 indicates the target STOP node was
%       not reached from the START node.  0 indicates the running
%       intersection property does NOT hold from START to STOP.
%
% See Cormen, Leiserson and Rivest, "An intro. to algorithms" 1994, p478.
%
% This file was adapted from matlabtools.googlecode.com
% by Michael McGeachie 2012.

n = length(tree);

global white gray black color
white = 0; gray = 1; black = 2;
color = white*ones(1,n);

global cycle breakflag runint
cycle = 0;
breakflag = false;
runint = true;

global pred
pred = zeros(1,n);

global sd
sd = false(1,n);

imems = unique([tree(start).members, tree(start).dmembers, tree(start).index],'legacy');
jmems = unique([tree(stop).members, tree(stop).dmembers, tree(stop).index],'legacy');
intset = intersect(imems, jmems, 'legacy');

dfs_visit(start, tree, stop, intset);

if (~breakflag)
  % searched to end of tree and didn't reach the STOP node
  rint = -1;
else
  if (runint)
    rint = 1;
  else
    rint = 0;
  end
end
hascycle = cycle;

end
%%%%%%%%%%

function dfs_visit(u, tree, stop, intset)

global white gray black color cycle pred sd breakflag runint
color(u) = gray;

% check intersection set:
umems = unique([tree(u).members, tree(u).dmembers, tree(u).index],'legacy');
sd(u) = isempty(setdiff(intset, umems,'legacy')); % must be empty for running set property to hold

if (u == stop)
  % set indicator to stop looking through children at higher depths.
  breakflag = true;
  if (~sd(u))
    runint = false;
  end
  color(u) = black;
  % stop looking at any further children
  return;
end

ns = cell2mat(tree(u).children);
ns = setdiff(ns, pred(u),'legacy'); % don't go back to visit the guy who called you!

for v=ns(:)'
  %fprintf('u=%d, v=%d, color(v)=%d\n', u, v, color(v))
  switch color(v)
    case white, % not visited v before
     pred(v)=u;
     dfs_visit(v, tree, stop, intset);
   case gray, % back edge - v has been visited, but is still open
    cycle = 1;
    %fprintf('cycle: back edge from v=%d to u=%d\n', v, u);
   case black, % v has been visited, but is closed
    % no-op
  end
  if (breakflag)
    if (~sd(u))
      runint = false;
    end
    break;
  end
end
color(u) = black;
end

