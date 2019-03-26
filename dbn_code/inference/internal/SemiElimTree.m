function tree = SemiElimTree(cs, eo)
% tree = SemiElimTree(cs, eo)
%
% Combines discrete nodes that are subsets of each other, taking a 
% Elimination Tree and then computing the SEMI ELIMINATION TREE from that tree.
% 
% These are both types of junction trees.  The Semi-Elimination tree is
% what we want for hybrid bayesian network inference
%
% INPUT : 
%   CS : array of ClusterSetTree objects
%   EO : Elimination ordering of NODES in the original bayes net
%
% OUTPUT :
%   TREE : an array of ClusterSetTree objects that forms a semi-elimination
%       tree, or junction tree of the bayesian network
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


tree = [];
list = cs;
i = 1;
while (i <= length(list))
    c = list(i);    
    % if cluster c is entirely discrete, check if it is a subset of other
    % discrete clusters:
    subset = false;
    if (list(i).discrete);
        for j = i+1:length(list)
            if (list(j).discrete)
                if (issubset(c, list(j)))
                    subset = true;
                    break;
                end
            end
        end 
    end        
    
    % if this is not a discrete subset of another discrete set, 
    % just add this to tree:
    if (~subset)
        if (isempty(tree))
            tree = c;
        else
            tree(end+1) = c;
        end;
        i = i + 1;
        continue;
    end
    
    % we do have a subset:
    % search these children in reverse elimination order,
    % so that the child chosen was eliminated after the other children
    bestchild = 0;
    minenum = 0;
    for j = 1:length(c.children)
        if (list(c.children{j}).discrete && issubset(c, list(c.children{j})))
            %enum = eo(list(c.children{j}).index);
            enum = find(eo == list(c.children{j}).index);
            if (enum > minenum)
                minenum = enum;
                bestchild = j;
            end
        end
    end
    
    % do stuff
    % this child will absorb c
    cp = list(c.children{bestchild});
    oldindex = c.children{bestchild};
    % set the parent:
    cp.parent = c.parent;
    % add children from the parent :
    for k = 1:length(c.children)
        if (k ~= bestchild)
            cp.children{end+1} = c.children{k};
        end
    end
    % set the parent of these children:
    for k = 1:length(cp.children)
        list(cp.children{k}).parent = i;
    end
    % this reorders one element of list while removing the i^th one:
    list = [list(1:i-1), cp, list(i+1:oldindex-1), list(oldindex+1:end)];
    % now update indices in the following list elements :
    for k = 1:length(list)
        if (list(k).parent >= oldindex)
            list(k).parent = list(k).parent-1;
        end
        for m = 1:length(list(k).children)
            if (list(k).children{m} >= oldindex)
                list(k).children{m} = list(k).children{m} -1;
            end
        end
    end
    % no need to increment the index i here, since we have instead
    % decreased the size of LIST by 1.
end


function subset = issubset(c, s)
% takes two CLUSTERSETS as input
n = [c.index, c.members];
d = [s.index, s.members];
%check subset:
subset = isempty(setdiff(n,d, 'legacy'));









