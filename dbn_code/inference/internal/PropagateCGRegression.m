function tree = PropagateCGRegression(tree, bn)
%tree = PropagateCGRegression(tree, bn)
% 
% Cowell's algorithm 5.2: initialization of conditional gaussian
% regressions in a hybrid bayesian network.
%
% INPUT : 
%   TREE : array of ClusterSetTrees representing the junction tree for the
%       bayesian network
%   BN : master array of NODES that is the original bayes net
%
% DEBUGGING STATE: Deal with cases where there are multiple exchanges 
% performed on the postbag.  However, this seems to never occur.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

useBioinfToolbox = true;

% make a top-sort of the NODES
if (useBioinfToolbox)
    % (this fucntion takes a sparse adjacency matrix)
    toporder = graphtopoorder(sparse(AdjMatrixNodes(bn)));
else
    toporder = TopSort(AdjMatrixNodes(bn));
end

for i = length(tree):-1:1
    
    % sort the "regressions" (aka LPPOTs) in the postbag
    % so their HEADs match the toporder earlier.
    for k = 1:length(tree(i).postbag)
        if (~isempty(tree(i).postbag{k}))
            tree(i).postbag{k} = sortpostbag(tree(i).postbag{k}, toporder);
        end
    end
    
    % cycle through each postbag and deal with each message
    for p = 1:length(tree(i).postbag)
        % now deal with all distributions in the postbag:
        ltp = length(tree(i).postbag{p});
        for k = 1:ltp;
            R = tree(i).postbag{p}(k);
            % if LPPot R has cst.index in the tail, then Exchange:
            if (ismember(tree(i).index, R.tail, 'legacy'))
                % lppotential is always length 1 -                
                if (isempty(tree(i).lppotential{p}))
                    % deal with this case:
                    error('Exchanging postbag with empty LP Potential list');
                end
                % do a sequence of EXCHANGES here for every member in the
                % postbag:
                if (length(tree(i).lppotential{p}) > 1)
                    % at this point, we should do the exchanges in
                    % Elimination Order of the nodes:
                    error('lp potential has too many elements');
                end
                newlp = tree(i).lppotential{p}(1);
                [newlp, R] = Exchange(newlp, R);
                % last remaining lppot is the one we retain in the node:
                tree(i).lppotential{p}(1) = newlp;
            end

            % put a copy of R in the postbag or lppotential of cst's parent:
            if (~isempty(tree(i).parent))
                if (~tree(tree(i).parent).discrete)
                    if (tree(tree(i).parent).index == R.head)
                        tree(tree(i).parent) = tree(tree(i).parent).AddToLPPotential(R);
                    else
                        tree(tree(i).parent) = tree(tree(i).parent).AddToPostbags(R);
                    end
                end
                % else : discrete nodes have been dealt with elsewhere.
            end
        end
        tree(i).postbag{p} = [];
    end
end


end

% sorts the contents of a postbag (LPPotentials) by some specific rank
% order
function sortbag = sortpostbag(pbag, rankorder)
    sortbag = pbag;
    found = false(size(pbag));
    for k = 1:length(sortbag)
        first = 0;
        for j = 1:length(pbag)
            if (~found(j))
                if (first == 0)
                    first = j;
                end
                if (rankorder(pbag(j).head) < rankorder(pbag(first).head))
                    first = j;
                end
            end
        end
        sortbag(k) = pbag(first);
        found(first) = true;
    end
end


