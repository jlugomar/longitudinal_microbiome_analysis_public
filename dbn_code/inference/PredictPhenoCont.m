function [acc, p, z, truelikelihood] = PredictPhenoCont(BN, verbose)
%[acc, p, z, truelikelihood] = PredictPhenoCont(BN, verbose)
%
% Main function for taking a given CG Bayesian Network that has had both
% its structure and its parameters learned, and performing inference:
% Predicting the phenotype of the data when the phenotype is a continuous
% variable.  Also works for discrete variables.  May take slightly longer
% than PredictPheno() for discrete variables, and that function could be
% used when likelihood computations are not required.
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
%   TRUELIKELIHOOD : log probability of observing a less likely outcome for
%       each variable than the value assigned to that var by the input
%       data.  
%
%
% Copyright Michael McGeachie, 2014.  MIT license. See cgbayesnets_license.txt.

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
lastdiscpred = [];
if (nodes(phenind).discrete)
    lastdiscpred = nodes(phenind).self;
end
% this is the list of column indices for all nodes, except the phenotype node
allnodes = false(size(colnames));
for i = 1:length(nodes)
    allnodes(nodes(i).colind) = true;
end
% this is a list of node indices for all nodes, except the phenotype:
evnodeinds = true(size(nodes));
evnodeinds(phenind) = false;
nlist = 1:length(nodes);
evnodeinds = nlist(evnodeinds);

% find the phenotype node in the tree:
for r = 1:length(tree)
    if (sum(strcmpi(tree(r).cpt.factors, nodes(phenind).self)) > 0)
        treeind = r;
        break;
    end
end


allnodes(phencol) = false;
ncases = size(data,1);
pred = -1 * ones(ncases, 1);
logpredprob = pred;
truelikelihood = pred;
p = ones(ncases,1) * -1;
z = p;
origtree = tree;
truedata = data(:,phencol);
tsd = std(truedata);
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
    lastpush = false;
    for j = 1:length(nodes)
        % don't enter evidence on the phenotype column
        if (j == phenind && nodes(j).discrete)
            tree = HybridDiscreteEvidence(tree, [], []);
        elseif (nodes(j).discrete)
            % this part just moves evidence from the continuous nodes to
            % the discrete nodes
            lastpush = true;
            tree = HybridDiscreteEvidence(tree, {nodes(j).self}, drow(nodes(j).colind));
        end
    end
    if (~lastpush)
        tree = HybridDiscreteEvidence(tree, {},[]);
    end
    % now do the discrete part
    [~, treeupdate] = DiscreteEvidence(tree, colnames(allnodes), ...
        drow(allnodes), lastdiscpred, evnodeinds, nodes);
    tree = treeupdate;

    % now do the prediction.
    if (nodes(phenind).discrete)    
        csttree = tree(treeind);
        % check the discrete cluster (the root of the junction tree) and it's CPT
        [maxval, mind] = max(csttree.cpt.logprob);
        logpredprob(i) = maxval;
        % report log probabilities:
        z(i) = exp(maxval - logplus(csttree.cpt.logprob));

        % find phenotype again : usually at last index:
        factorind = find(strcmp(csttree.cpt.factors,nodes(prednodeinds(j)).self));
        p(i) = csttree.cpt.values{factorind}(mind);
        if (verbose)
            fprintf(1,'\t prediction for example %d is : %d with confidence %1.2f ',i,p(i), z(i));
            if (drow(phencol) == p(i))
                fprintf(1,' right!\n');
            else
                fprintf(1,' wrong!\n');
            end
        end
    else
        % compute expected value of continuous phenotype node:
        % also checks the prob of the actual data:
        [p(i), logpredprob(i), truelikelihood(i)] = ...
            RecursiveContExp(tree,phenind,[],drow(nodes(phenind).colind));
        z(i) = logpredprob(i);
        if (verbose)
            fprintf(1,'\t prediction for example %d is : %2.2f . True value is : %2.2f',i, p(i), drow(phencol));
            
            % report error in avg z-score:
            zmiss = abs(p(i) - drow(phencol))./tsd;
            fprintf(' . Z-error : %2.2f\n',zmiss);
        end
        
    end
    
end

% compute final average performance metrics:
if (nodes(phenind).discrete)
    acc = (p == data(:,phencol));
    acc = sum(acc) / ncases;
else
    % report error in z-score:
    acc = abs(p - truedata)./tsd;
end

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
