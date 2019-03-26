function [acc, p, z] = PredictPheno(BN, verbose)
%[acc, p, z] = PredictPheno(BN, verbose)
%
% Main function for taking a given CG Bayesian Network that has had both
% its structure and its parameters learned, and performing inference:
% Predicting the phenotype of the data.
%
%  INPUT: 
%  BN : a BayesNet class object that represents the Bayes Net to predict
%       upon.  Must have these fields appropriate instantiated:
%    BN.pheno : string indicating which column to predict
%    BN.data : a matrix of data to predict
%    BN.cols : cell array of strings for column titles, or variable names.
%       must match BN.data and one element matches BN.pheno
%    BN.nodes : master list of NODES making up the bayes net, instantiated
%       by LearnParams()
%    BN.tree : ClusterSetTree instantiated by LearnParams()
%  VERBOSE : if true increases output. default = false;
%
%  OUTPUT:
%   acc : the accuracy in % correct predictions of examples in FILENAME
%   p : the predicted phenotype values for each example in FILENAME
%   z : the prediction confidence of the predictions p.  between (0,1); but
%       normalized so that binary predictions are always in (0.5, 1).
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    verbose = false;
end

colnames = BN.cols;
nodes = BN.nodes;
tree = BN.tree;
pheno = BN.pheno;
data = BN.data;

% need to double check that the column names map appropriately to the nodes
%nodes = AddColIndstoNodes(nodes, colnames);

% phenotype index:
[phenind, phencol] = GetPhenoInds(pheno, nodes, colnames);
% phenind is the index for NODES(PHENIND).self == PHENO
% phencol is the index for DATA(:,PHENCOL)

% this is the list of indices for all nodes, except the phenotype node
allnodes = false(size(colnames));
for i = 1:length(nodes)
    allnodes(nodes(i).colind) = true;
end
allnodes(phencol) = false;

found = false;
% find the phenotype node in the tree:
for r = 1:length(tree)
    for rt = 1:length(tree(r).cpt.factors)
        if (sum(strcmpi(tree(r).cpt.factors{rt}, nodes(phenind).self)) > 0)
            treeind = r;
            found = true;
            break;
        end
        if (found)
            break;
        end
    end
%    if (sum(strcmpi(tree(r).cpt.factors, nodes(phenind).self)) > 0)
%        treeind = r;
%        break;
%    end
end

ncases = size(data,1);
p = ones(ncases,1) * -1;
z = p;
origtree = tree;
%parfor i = 1:ncases
for i = 1:ncases
     tree = origtree;
    % this assignment makes it clear to MATLAB's paralleltoolkit that DATA
    % can be split by row and sent to each process separately:
    drow = data(i,:);
    for j = 1:length(nodes)
        % don't enter evidence on the phenotype column
        if (j == phenind)
            continue;
        end
        if (~nodes(j).discrete)
            % read continuous values from datafile:
            tree = PushEvidence(tree, j, drow(nodes(j).colind));
        end
    end
    for j = 1:length(nodes)
        % don't enter evidence on the phenotype column
        if (j == phenind && nodes(j).discrete)
            tree = HybridDiscreteEvidence(tree, [], []);
        elseif (nodes(j).discrete)
            % this part just moves evidence from the continuous nodes to
            % the discrete nodes
            tree = HybridDiscreteEvidence(tree, {nodes(j).self}, drow(nodes(j).colind));
        end
    end
    % now do the discrete part
    [~, tree] = DiscreteEvidence(tree, colnames(allnodes), ...
        drow(allnodes), nodes(phenind).self);

    csttree = tree(treeind);

    % now do the prediction.
    % check the phenotype cluster (the root of the junction tree) and it's CPT
    [maxval, mind] = max(csttree.cpt.logprob);
    z(i) = exp(maxval - logplus(csttree.cpt.logprob));
    % predcition is the highest value:
    % find phenotype again : usually at last index:
    if (strcmpi(csttree.cpt.factors{length(csttree.cpt.factors)}, nodes(phenind).self))
        rootphenoind = length(csttree.cpt.factors);
    else
        % can be an error for continuous networks;
        %error('Conditional Probability Table has reorganized factors');
        
        % for discrete networks just find the index:
        for r = 1:length(csttree.cpt.factors)
            if (strcmpi(csttree.cpt.factors{r}, nodes(phenind).self))
                rootphenoind = r;
                break;
            end
        end
    end
    
    p(i) = csttree.cpt.values{rootphenoind}(mind);
    if (verbose)
        fprintf(1,'\t prediction for example %d is : %d with confidence %1.2f ',i,p(i), z(i));
        if (drow(phencol) == p(i))
            fprintf(1,' right!\n');
        else
            fprintf(1,' wrong!\n');
        end
    end
end

acc = (p == data(:,phencol));
acc = sum(acc) / ncases;

end

%% helper function for picking out phenotype indices
function [phenind, phencol] = GetPhenoInds(pheno, nodes, colnames)
phenind = 0;
phencol = 0;
if (ischar(pheno))
    for i = 1:length(nodes)
        if (strcmpi(nodes(i).self, pheno))
            phenind = i;
        end
        for j = 1:length(colnames)
            if (strcmpi(colnames{j}, pheno))
                phencol = j;
                break;
            end
        end
    end
else
    phencol = pheno;
    for i = 1:length(nodes)
        if (pheno == nodes(i).colind)
            phenind = i;
            break;
        end
    end
end
end
