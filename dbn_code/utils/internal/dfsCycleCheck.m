function hascycle = dfsCycleCheck(adjmat, start)
% DFS Perform a depth-first search of the graph starting from 'start'.
% hascycle = dfsCycleCheck(adj_mat, start)
%
% Input:
% 'adjmat' is an adjacency matrix of the network
% 'start' is the index into network to start the search 
%
% Output:
% 'hascycle' is true iff a directed cycle is found.
%
% See Cormen, Leiserson and Rivest, "An intro. to algorithms" 1994, p478.
%
% This file was adapted from matlabtools.googlecode.com
% by Michael McGeachie 2012.

n = length(adjmat);

global white gray black color
white = 0; gray = 1; black = 2;
color = white*ones(1,n);

global cycle breakflag
cycle = 0;
breakflag = false;

global pred
pred = zeros(1,n);

dfs_visit(start, adjmat);
hascycle = cycle;

end
%%%%%%%%%%

function dfs_visit(u, adjmat)

global white gray black color cycle pred breakflag
color(u) = gray;

if (breakflag)
  return;
end

% list of all children 
children = find(adjmat(u,:) > 0);

for v=children(:)'
  %fprintf('u=%d, v=%d, color(v)=%d\n', u, v, color(v))
  switch color(v)
    case white, % not visited v before
     pred(v)=u;
     dfs_visit(v, adjmat);
   case gray, % back edge - v has been visited, but is still open
    cycle = 1;
    breakflag = true;
    %fprintf('cycle: back edge from v=%d to u=%d\n', v, u);
   case black, % v has been visited, but is closed
    %no-op
  end
  if (breakflag)
    break;
  end
end
color(u) = black;
end

