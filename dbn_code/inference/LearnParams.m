function [treeOut, nodes] = LearnParams(nodesketch, filename, data, colnames, ...
    priorPrecision, discvals, is_discrete)
%[treeOut, nodes] = LearnParams(nodesketch, filename, data, colnames, priorPrecision, discvals, is_discrete)
%
% Main function for learning the parameters of a CG Bayesian Network. Use
% after learning the structure of the Bayesian Network.
%
%  INPUT: 
%  NODESKETCH: array of linked NODE class instantiations :
%  FILENAME : string pointing to the datafile to read and learn from. Can
%   be an empty array; will then use the DATA and COLNAMES input next.
%  DATA: matrix of data observations; can be used instead of a filename.
%  COLNAMES: list of column names; can be used instead of a filename
%  PRIORPRECISION: structure containing information about the bayesian
%   priors to use in parameter estimation
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%  DISCVALS : optional parameter; cell array of cells for values for every
%       discrete variable.  Use when there may be variable values
%       unrepresented in the input DATA.
%  IS_DISCRETE : optional parameter, an array of true/false indicating if
%       each node is discrete or not
%
%  OUTPUT:
%  TREEOUT: array of ClusterSetTree class objects representing the CG Bayes Net.
%  NODES: array of NODE class objects that represent nodes in the CG Bayes Net.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

% for testing and timing:
VERBOSE = false;

if (nargin < 3)
    [data,colnames] = ReadHeaderDataFile(filename);
end
if (nargin < 5)
    priorPrecision.nu = 1;
    priorPrecision.alpha = 1;
    priorPrecision.sigma = 1;
end
if (VERBOSE)
    fprintf(1,'computing NODETREE...\n');
end
if (nargin >= 7)
    bn = nodetree(nodesketch, data, colnames, discvals, is_discrete);
elseif (nargin >= 6)
    bn = nodetree(nodesketch, data, colnames, discvals);
else
    bn = nodetree(nodesketch, data, colnames);
end
if (VERBOSE)
    fprintf(1,'computing ELIMORDERING...\n');
end
eo = ElimOrdering(bn);
if (VERBOSE)
    fprintf(1,'computing MORALIZE...\n');
end
nodes = Moralize(bn);
if (VERBOSE)
    fprintf(1,'computing TRIANGULATE...\n');
end
nodes = Triangulate(nodes, eo);
if (VERBOSE)
    fprintf(1,'computing FORMCLUSTERSETS...\n');
end
cs = FormClusterSets(nodes, eo);
if (VERBOSE)
    fprintf(1,'computing ELIMTREE...\n');
end
cs = ElimTree(cs, eo);
if (VERBOSE)
    fprintf(1,'computing SEMIELIMTREE...\n');
end
tree = SemiElimTree(cs, eo);
%CheckRunningIntersection(tree);
if (VERBOSE)
    fprintf(1,'computing INITCLUSTERSETTREE...\n');
end
tree = InitClusterSetTree(tree, nodes);
if (VERBOSE)
    fprintf(1,'computing INITCLUSTERPOTENTIAL...\n');
end
tree = InitClusterPotential(tree, nodes, data, priorPrecision);
if (VERBOSE)
    fprintf(1,'computing PROPAGATECGREGRESSION...\n');
end
treeOut = PropagateCGRegression(tree, bn);


