function [adjmatNew, weightMatrix] = AddNodeToBN(adjmat, data, cols, pheno, addNode, ...
    priorPrecision, onlyMB, weightMatrix)
%adjmatNew = AddNodeToBN(adjmat, data, cols, pheno, addNode, priorPrecision, onlyMB)
%
% add a link from ADDNODE to the existing bayesian network in ADJMAT.
% Searches for best link from ADDNODE to existing network, which possibly
% is not in the Markov Blanket of the phenotype node.
%
% INPUT:
% ADJMAT: adjacency matrix representing the current Bayes Net.
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% ADDNODE: a string representing the node to be added to the Bayes Network.
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% ONLYMB: if true, will only add the ADDNODE in a way that it contributes
%   to the markov blanket of the phenotype node.  Default = true.
% WEIGHTMATRIX: a matrix of edge weights for each edge added to the
%   network.  optional.
%
% OUTPUT:
% ADJMATNEW: new adjacency matrix that represents the network with the
%   ADDNODE added to the network.
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

if (nargin < 7)
    onlyMB = true;
end
if (nargin < 8)
    weightMatrix = adjmat;
end

hinds = sum(adjmat,1) > 0;
vinds = sum(adjmat,2) > 0;
inds = hinds | vinds';

f = strcmpi(pheno, cols);
phnind = find(f);
f = strcmpi(addNode, cols);
addind = find(f);

disc = IsDiscrete(data);

% set up a matrix of possible edges to add:
possible = zeros(size(adjmat));
if (onlyMB)
    possible(phnind, addind) = 1;
    % find all children of phenotype:
    children = adjmat(phnind, :) == 1;
    possible(addind,children) = ones(1,sum(children));
else
    possible(inds, addind) = ones(sum(inds),1);
    possible(addind, inds) = ones(1,sum(inds));
    if (disc(addind))
        % can't have links from non disc to this
        possible(~disc,addind) = zeros(sum(~disc),1);
    else
        % can't have links from this to disc
        possible(addind, disc) = zeros(1,sum(disc));
    end
    % can't have links to the phenotype from anything:
    possible(addind, phnind) = 0;
end

lldiff = -Inf * ones(length(disc));
for i = 1:length(disc)
    for j = 1:length(disc)
        if (possible(i,j))
            %% consider making an edge from node i to j.

            % data of the child node
            childData = data(:,j);

            % initial model of node j 
            parents = adjmat(:,j) ~= 0;
            contParentData = data(:,parents' & ~disc);
            discParentData = data(:,parents' & disc);
            if (disc(j))
                % is a discrete node
                model_ini=learnLocalBN_DiscToDisc(discParentData',childData',priorPrecision);
            else
                model_ini=learnLocalBN_MixToCont(contParentData',discParentData',childData',priorPrecision);
            end
            % now add edge (i,j)
            parents2 = parents;
            parents2(i) = true;
            contParentData = data(:,parents2' & ~disc);
            discParentData = data(:,parents2' & disc);
            if (disc(j))
                % is a discrete node
                model_addedge=learnLocalBN_DiscToDisc(discParentData',childData',priorPrecision);
            else
                model_addedge=learnLocalBN_MixToCont(contParentData',discParentData',childData',priorPrecision);
            end

            lldiff(i,j) = model_addedge.logLLH - model_ini.logLLH;
        end
    end
end

% Find the best edge:
[~,maxi] = max(max(lldiff,[],2));
[bestdiff,maxj] = max(lldiff(maxi,:));

% go with it
adjmatNew = adjmat;
adjmatNew(maxi,maxj) = 1;

% and keep the weight : 
weightMatrix(maxi,maxj) = bestdiff;



