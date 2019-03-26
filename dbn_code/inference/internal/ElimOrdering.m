function e = ElimOrdering(nodes, usetopsort)
%e = ElimOrdering(nodes, usetopsort)
%
% generate an ELIMINATION ORDERING from a list of nodes, where nodes is an
% array of NODE structs with:
% node(n).discrete: TRUE if node is discrete, FALSE if node is continuous
%
% elimination ordering has all continuous vars before all discrete vars.
% uses a reverse topological sort on the continuous nodes, and concatenates
% it with a reverse topological sort on the discrete nodes
%
% TOOLBOXES: can set the constant useBioinfToolbox = false, and it will not
% use the MATLAB bioinformatics toolbox.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.
%

if (nargin < 2)
    usetopsort = true;
end

useBioinfToolbox = true;

discnodes = [];
contnodes = [];

for i = 1:length(nodes)
    if (nodes(i).discrete)
        discnodes(end+1) = i;
    else
        contnodes(end+1) = i;
    end
end


if (usetopsort)
    % get adjacency matrices:
    contadj = AdjMatrixNodes(nodes, false);
    discadj = AdjMatrixNodes(nodes, true);

    % then call graphtopoorder
    % note that these node indices must be in sorted order

    if (useBioinfToolbox)
        di = graphtopoorder(sparse(discadj(discnodes, discnodes)));
        ci = graphtopoorder(sparse(contadj(contnodes, contnodes)));
        e = [contnodes(ci(end:-1:1)), discnodes(di(end:-1:1))];
    else
        di = TopSort(discadj(discnodes, discnodes));
        ci = TopSort(contadj(contnodes, contnodes));
        e = [discnodes(di), contnodes(ci)];
        e = e(end:-1:1);
    end
    
else
    e = [contnodes, discnodes];
end



