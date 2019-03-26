function [fprod, tree] = DiscreteEvidence(tree, allvarnames, evidence, prodvar, evnodeinds, nodes)
%[fprod, tree] = DiscreteEvidence(tree, allvarnames, evidence, prodvar, evnodeinds, nodes)
%
% function for entering discrete nodes evidence to the discrete part of a
% hybrid discrete/continuous Bayesian network.
%
% INPUT : 
%   TREE : the junction tree, an array of ClusterSetTree objects
%   ALLVARNAMES : cell array of names of variables for which evidence has
%       been observed.
%   EVIDENCE : the values of ALLVARNAMES that have been observed
%   PRODVAR : the product variable name, a string representing which node
%       to perform the factor product operation upon.
%   EVNODEINDS : indices into NODES that indicate discrete nodes for which
%       evidence has been observed.  These are propagated onto the continuous
%       nodes. Optional.  If empty, no discrete evidence is propagated to
%       continuous nodes.
%   NODES : the cannonical NODE array for the BN; required if EVNODEINDS is
%       used.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    evnodeinds = [];
end
if (nargin < 6)
    nodes = [];
end

% collect all CPTs from each discrete ClusterSetTree
cpts = repmat(CondProbTable(), size(tree));
cptindex = false(size(tree));
for i = 1:length(tree)
    if (tree(i).discrete)
        cpts(i) = tree(i).cpt;
        cptindex(i) = true;
    end
end


% step 2: reduce all conditional probability tables, to 
% their evidence
freds = factorreduce(cpts(cptindex),allvarnames,evidence);

% step 3: compute a product of these factors on PRODVAR
fprod = factorproduct(freds, prodvar);
if (isempty(fprod))
    fprod = freds;
end

% put these back into the tree, where we found them:
factorind = 1;
if (length(fprod) == sum(cptindex))
    for i = 1:length(cptindex)
        if (cptindex(i))
            tree(i).cpt = fprod(factorind);
            factorind = factorind + 1;
        end
    end
else
    for i = 1:length(cptindex)
        if (cptindex(i))
            tree(i).cpt = fprod(1);
        end
    end
end

% step 4: enter discrete evidence into all continuous clustersets:
if (~isempty(evnodeinds))
    disclist = false(size(evnodeinds));
    evindex = zeros(size(evnodeinds));
    for i = 1:length(evnodeinds)
        evindex(i) = find(strcmp(nodes(evnodeinds(i)).self,allvarnames));
        if (nodes(evnodeinds(i)).discrete)
            disclist(i) = true;
        end
    end
    for i = 1:length(tree)
        if (~tree(i).discrete)
            % this uses node indices, rather than varnames.
            tree(i) = tree(i).DiscConditionEvidence(evnodeinds(disclist), ...
                evidence(evindex(disclist)));
        end    
    end
end
