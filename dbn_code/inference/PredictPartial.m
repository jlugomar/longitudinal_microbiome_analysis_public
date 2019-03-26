function [pred,z,logevprob,truelikelihood] = PredictPartial(BN, evvars, predvars)
%[pred,z,logevprob,truelikelihood] = PredictPartial(BN, evvars, predvars)
%
% General function for using a BayesNet to predict some of the variables
% given the evidence on the other variables.  
%
% Also predicts continuous phenotypes.
%
%  INPUT: 
%  BN : a BayesNet class object that represents the Bayes Net to predict
%       upon.  Must have these fields appropriate instantiated:
%    BN.data : a matrix of data to predict
%    BN.cols : cell array of strings for column titles, or variable names.
%       must match BN.data and one element matches BN.pheno
%    BN.nodes : master list of NODES making up the bayes net, instantiated
%       by LearnParams()
%    BN.tree : ClusterSetTree instantiated by LearnParams()
%  EVVARS : Variables for which evidence has been observed.  Values of
%       BN.data(:,EVVARS) will be used.
%  PREDVARS : Variables that will be predicted.  Input as indices into
%       BN.cols.
%
%  OUTPUT:
%   PREDS : the predicted values for each PREDVARS given each row of
%       evidence assigned to the EVVARS. Of size [length(BN.data(:,1)) by
%       length(PREDVARS)]
%   Z : the prediction confidence of the predictions p.  between (0,1); but
%       normalized so that binary predictions are always in (0.5, 1).  Not
%       defined for predictions of continuous variables.
%   LOGEVPROB : log probability of the evidence; only defined for
%       predictions of entirely discrete variables.
%   TRUELIKELIHOOD : log probability of observing a less likely outcome for
%       each variable than the value assigned to that var by the input
%       data.  
%
%
% Copyright Michael McGeachie, 2014.  MIT license. See cgbayesnets_license.txt.

colnames = BN.cols;
nodes = BN.nodes;
tree = BN.tree;
data = BN.data;

% evvars are inds and so are predvars, need to convert to node inds:
evnodeinds = zeros(size(evvars));
prednodeinds = zeros(size(predvars));
for i = 1:length(evvars)
    evnodeinds(i) = GetPhenoInds(evvars(i),nodes, colnames);
end
lastdiscpred = [];
for i = 1:length(predvars)
    prednodeinds(i) = GetPhenoInds(predvars(i),nodes, colnames);
    if (nodes(prednodeinds(i)).discrete)
        lastdiscpred = nodes(prednodeinds(i)).self;
    end
end

ncases = size(data,1);
origtree = tree;
pred = -1 * ones(ncases, length(predvars));
logpredprob = pred;
truelikelihood = pred;
z = pred;
logevprob = -1 * ones(ncases,1);
%parfor i = 1:ncases
for i = 1:ncases
    tree = origtree;
    % this assignment makes it clear to MATLAB's paralleltoolkit that DATA
    % can be split by row and sent to each process separately:
    drow = data(i,:);
    for j = 1:length(evnodeinds)
        % enter all evidence from the nodes we have evidence on, starting
        % with continuous ones:
        if (~nodes(evnodeinds(j)).discrete)
            tree = PushEvidence(tree, evnodeinds(j), drow(nodes(evnodeinds(j)).colind));
        end
    end
    lastpush = false;
    for j = 1:length(evnodeinds)
        % now do the discrete evidence:
        if (nodes(evnodeinds(j)).discrete)
            % this part just moves evidence from the continuous nodes to
            % the discrete nodes
            lastpush = true;
            tree = HybridDiscreteEvidence(tree, {nodes(evnodeinds(j)).self}, drow(nodes(evnodeinds(j)).colind));
        end
    end
    % make sure we call HybridDiscreteEvidence even if we don't have any
    % discrete evidence to enter:
    if (~lastpush)
        tree = HybridDiscreteEvidence(tree, {},[]);
    end
    % now do the discrete part
    [~, treeupdate] = DiscreteEvidence(tree, colnames(evvars), ...
        drow(evvars), lastdiscpred, evnodeinds, nodes);
    tree = treeupdate;
    
    % can read the evidence likelihood from the normalization constant in
    % the head CST:
    if (isempty(tree(1).cpt)) 
        % no discrete nodes results in this kind of cpt:
        logevprob(i) = logplus(tree(1).logweighttable);
    else
        % this isn't the probability... it's an unnormalized likelihood (max of the PDF):
        % logevprob(i) = logplus(tree(1).cpt.logprob);
        logevprob(i) = 0;
    end
    
    %% predict with only the nodes requested for prediction.
    for j = 1:length(prednodeinds)
        % predict for each prediction node:
        % find the factor associated with that node:
        % for discrete networks just find the index:
        if (nodes(prednodeinds(j)).discrete)
            for r = 1:length(tree)
                if (sum(strcmpi(tree(r).cpt.factors, nodes(prednodeinds(j)).self)) > 0)
                    treeind = r;
                    break;
                end
            end
            csttree = tree(treeind);
            % check the discrete cluster (the root of the junction tree) and it's CPT
            [maxval, mind] = max(csttree.cpt.logprob);
            logpredprob(i,j) = maxval;
            % report (actual, non-log) probabilities:
            z(i,j) = exp(maxval - logplus(csttree.cpt.logprob));
            % need to find the index within this factor for the value:
            % look through csttree.cpt.factors, find index of one that
            % matches what we're predicting: 
            factorind = find(strcmp(csttree.cpt.factors,nodes(prednodeinds(j)).self));
            pred(i,j) = csttree.cpt.values{factorind}(mind);
        else
            % also checks the prob of the actual data:
            [pred(i,j), logpredprob(i,j), truelikelihood(i,j)] = ...
                RecursiveContExp(tree,prednodeinds(j),[],drow(nodes(prednodeinds(j)).colind));
            z(i,j) = logpredprob(i,j);

        end
    end
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
