function [expct, probextreme, hprobnorm] = RecursiveContExp(tree, var, condind, hypoval)
% [expct, probextreme, hprobnorm] = RecursiveContExp(tree, var, condind, hypoval)
% 
% this computes the expected value of a continuous Gaussian node given a
% ClusterSetTree and the conditioning index of any discrete nodes this
% continuous node may depend upon.  Recursively calls itself to evaluate
% conditional expected values of continuous nodes that it depends upon.
%
%
% INPUT :
%   TREE : ClusterSetTree array representing the Bayes Net
%   VAR : the variable to predict, an index into the canonical NODE array
%       for the BN.
%   CONDIND : the indices of conditioning discrete node values, into 
%       TREE.lppotential{CONDIND} and TREE.cpt.logprob(CONDIND). Can be empty.  
%       If empty, computes expectation across all possible conditioning 
%       assignments to discrete parents according to the probabilities of
%       those assignments.
%   HYPOVAL : hypothetical value for NODE(VAR), will output the probability of
%       seeing a value more extreme than HYPOVAL for VAR in HPROBNORM
%
% OUTPUT :
%   EXPCT : the expected value of NODE(VAR), given discrete vars are set to
%   PROBEXTREME : the probability of seeing a value more extreme that the
%       one predicted for NODE(VAR).  Usually = 0.5 since we predict expected
%       value.
%   HPROBNORM : the log probablity of seeing a value more extreme than
%       HYPOVAL for NODE(VAR)
%
%
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.



if (nargin < 3)
    condind = [];
end
if (nargin < 4)
    hypoval = [];
end

for r = 1:length(tree)
    if (tree(r).index == var)
        treeind = r;
        break;
    end
end

% TODO: make sure that the one discrete cluster in the tree is tree(1)
% or, if it's not discrete just don't care about it.
nodisc = false;
if (~tree(1).discrete)
    nodisc = true;
end
% always just 1 discrete tree cluster?
%for i = 1:length(tree)
%    if (i > 1 && tree(i).discrete)
%        error('Expecting exactly one discrete CST');
%    end
%end

% find the discparenttree:
i = r;
while (~tree(i).discrete)
    if (~isempty(tree(i).parent))
        i = tree(i).parent;
    else
        nodisc = true;
        break;
    end
end
discparenttree = i;

csttree = tree(treeind);
v = zeros(size(csttree.lppotential));
w = v;
cdfprob = v;
hypoprob = v;
if (isempty(condind))
    range = 1:length(csttree.lppotential);
else
    range = condind;
end
for k = range
    % if this regression is conditioned on some other vars, we
    % have to compute THOSE expectations, since by linearity of
    % expectatino we have E[v] = const + sum_i(param(i) * E[tail(i)])
    klpp = csttree.lppotential{k};
    paramval = zeros(size(klpp.params));
    for i = 1:length(csttree.lppotential{k}.params)
        paramval(i) = RecursiveContExp(tree, csttree.lppotential{k}.tail(i), k);
        klpp = klpp.AddEvidence(csttree.lppotential{k}.tail(i),paramval(i));
    end
    % should now be an unconditional regression:
    v(k) = klpp.const;
    % uses log probs :
    [w(k),cdfprob(k)] = klpp.MakeWeight(klpp.head, v(k));
    % need to see if values of discrete vars have been entered as evidence,
    % in which case, we need to limit our attention to just those cases.
    
    % could be no discrete node:
    if (~nodisc)
        w(k) = w(k) + tree(discparenttree).cpt.logprob(k);
    end
    if (~isempty(hypoval))
        [~,hypoprob(k)] = klpp.MakeWeight(klpp.head, hypoval);
    end
end
% normalize probabilities here:
totweights = logplus(w(range));
wnormed = w(range) - totweights;
% have to do this in non-log prob since we multiple it by values to find
% expected value:
wnormed = exp(wnormed);
% compute expectation E(prednodeinds(j)) :


% normalize cdfprob() here:
if (~nodisc) 
    % if there's no discrete conditions, don't have to take expectation across cont vars:
    % these aren't true probabilities, these are relative likelihoods,
    % which are in general > 1.  They are then normalized to be
    % probabilities:
    
    % instead of using tree(1) here, should probably follow the csttree up
    % the chain until we find the first discrete parent CST:
    
    elogprob = tree(discparenttree).cpt.logprob(range);
    enormprob = logplus(tree(discparenttree).cpt.logprob(range));
    cdfprobnorm = logplus(cdfprob(range) + (elogprob - enormprob));
else
    % this is length 1
    if (length(cdfprob) ~= 1)
        error('Error: no discrete conditions, but found multiple conditioning vars');
    end
    cdfprobnorm = cdfprob;
end

% normalize hypoprob() here :
hprobnorm = -1;
if (~isempty(hypoval))
    if (~nodisc)
        hprobnorm = logplus(hypoprob(range)+ (elogprob - enormprob));
    else
        hprobnorm = hypoprob;
    end
end

probextreme = cdfprobnorm;
expct = v(range)' * wnormed;












