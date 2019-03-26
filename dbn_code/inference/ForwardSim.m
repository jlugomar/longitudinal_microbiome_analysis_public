function [predictions, likelihood, truelikelihood] = ...
    ForwardSim(BN, CSTree, numrecurse, recursiveInputCallback)
% [predictions, likelihood, truelikelihood] = ForwardSim(BN, CSTree, numrecurse, recursiveInputCallback)
% 
% take a BayesNet object and follow the edges forward, rather than
% backwawrd (as in typical inference) to arrive at the most likely values
% for each variable, given the values of its parents.  This is useful for
% Dynamic Bayes Nets that only have edges going from prior time points to
% future time points.
%
% INPUT :
%   BN : BayesNet object
%   CSTSREE : A ClusterSetTree array, generated with output from
%       ForwardSimLearnParams().  Optional.  If not provided, will call
%       ForwardSimLearnParams().
%   NUMRECURSE : number of times to do recursive input of prediction.
%       Starts with BN.data(1,:) and then fills up to
%       predictions(1:numrecurse,:)
%   RECURSIVEINPUTCALLBACK : function pointer to a function that
%       post-processes prediction data and turns it into input data.
%       Necessary for recursive prediction if, eg, normalization is
%       required.  Syntax: newinput = @recursiveInputCallback(inputarray);
%
% OUTPUT :
%   PREDICTIONS : values for each node in the BN
%   LIKELIHOOD : log prob of predicted values
%   TRUELIKELIHOOD : log prob of actual values
%
%
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

if (nargin < 4)
    recursiveInputCallback = [];
end
if (nargin < 3)
    numrecurse = 0;
end
if (nargin < 2)
    cs = ForwardSimLearnParams(BN);
else
    cs = CSTree;
end

% to do :
% limit prediction to a certain subset of variables -

nodes = BN.nodes;

origcs = cs;
predictions = zeros(length(BN.data(:,1)),length(cs));
likelihood = -1 * ones(length(BN.data(:,1)),length(cs));
truelikelihood = -1 * ones(length(BN.data(:,1)),length(cs));
if (numrecurse == 0)
    limit = length(BN.data(:,1));
    input = BN.data;
else
    limit = numrecurse;
    input = BN.data(1,:);
end
nodetocstmap = zeros(size(nodes));
for i = 1:length(nodes)
    for j = 1:length(cs)
        if (cs(j).index == i)
            nodetocstmap(i) = j;
            break;
        end
    end
end
for d = 1:limit
    cs = origcs;
    % predict on each node, enter evidence on all neighbors of that node:
    for i = 1:length(nodes)
        n = nodes(i);
        % just read the neighbors from the cs(j):
        neighbors = cs(nodetocstmap(i)).members;
        % make sure we remove itself from this list:
        % (dmembers will include the self node, but members will not?)
        neighbors = setdiff(neighbors, nodetocstmap(i));
        % this is similar to the sequence in PushEvidence():
        if (~n.discrete)
            % continuous node
            for k = 1:length(neighbors)
                % do all the continuous evidence
                if (~nodes(nodetocstmap(neighbors(k))).discrete)
                    var = nodetocstmap(neighbors(k));
                    varval = input(d,nodes(nodetocstmap(neighbors(k))).colind);
                    for j = 1:length(cs(nodetocstmap(i)).lppotential)
                        cs(nodetocstmap(i)).lppotential{j}(1) = ...
                            cs(nodetocstmap(i)).lppotential{j}(1).AddEvidence(var, varval);
                    end
                end
            end
            
            % then do the discrete evidence:            
            dmems = cs(nodetocstmap(i)).dmembers;
            if (~isempty(dmems))
                % discrete evidence is just selecting the right
                % lppotential based on the conditioning variables'
                % values

                % 1) all the discrete evidence simultaneously; 
                vars = dmems;
                varcols = vars;
                for k = 1:length(vars)
                    varcols(k) = nodes(vars(k)).colind;
                end
                varassigns = input(d,varcols);
                valinds = zeros(size(varassigns));
                for v = 1:length(valinds)
                    for vi = 1:length(cs(nodetocstmap(i)).dvalues{v})
                        if (varassigns(v) == cs(nodetocstmap(i)).dvalues{v}{vi})
                            valinds(v) = vi;
                            break;
                        end
                    end
                end
                ind = cs(nodetocstmap(i)).dvalstepsize * (valinds -1)' + 1;
                % 2) just use that to index into lppotential{} and
                if (ind <= 0)
                    warning('index <= 0, index = %d',ind);
                end
            else
                ind = 1;
            end
            % 3) read out the lppotential.const as the value prediction.
            predictions(d,n.colind) = cs(nodetocstmap(i)).lppotential{ind}.const;
            % also check the prob of the actual data:
            [~,truelikelihood(d,n.colind)] = ...
                cs(nodetocstmap(i)).lppotential{ind}.MakeWeight(n.index, input(d,n.colind));
        else
            % do discrete prediction based on parent assignments -
            dmems = cs(nodetocstmap(i)).dmembers;
            % trick thing here: node includes itself as a discrete member,
            % but we don't want to condition on the current value of
            % itself.  
            vars = setdiff(dmems, cs(nodetocstmap(i)).index);
            varcols = vars;
            for k = 1:length(vars)
                varcols(k) = nodes(vars(k)).colind;
            end
            
            % just use factorreduceworker(), which does what we want.  Except that
            % we need to give it variable NAMES (strings) for the vars.
            varnames = cell(1,length(vars));
            for v = 1:length(vars)
                varnames{v} = cs(nodetocstmap(vars(v))).name;
            end
            varassigns = input(d,varcols);
            [newcpt, ~] = factorreduceworker(cs(nodetocstmap(i)).cpt, varnames, varassigns, cs(nodetocstmap(i)));
            
            % just take max over : (cs(i).cpt.logprob), get index from that
            [likelihood(d,n.colind),mind] = max(newcpt.logprob);
            predictions(d,n.colind) = newcpt.values{1}(mind);
            truevalind = newcpt.values{1} == input(d,cs(nodetocstmap(i)).index);
            [truelikelihood(d,n.colind)] = newcpt.logprob(truevalind);
        end
    end
    if (numrecurse > 0)
        % should use a callback here to allow custom normalization of input
        % streams
        inline = recursiveInputCallback(predictions(d,:));
        input = [input;inline];
    end
end







