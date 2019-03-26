function cst = HybridDiscreteEvidence(cst, varnames, evvals)
%cst = HybridDiscreteEvidence(cst, varnames, evvals)
%
%  this function pushes evidence from the continuous part of a hybrid bayesian
%  network up into the discrete part of that network.
%
% Input:
%   CST : Array of ClusterSetTree class objects, describing the junction
%       tree.
%   VARNAMES : variable for which evidence is entered
%   EVVALS : value for variable VARNAME.
%
% DEBUGGING : Does this work on long chains of discrete ClusterSets?.
% ALso: not yet implemented for more than one variable at a time, which is
% probably fine.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% find boundary nodes, starting from the end:
for n = length(cst):-1:2

    % operate upon nodes that have non-empty weight tables, 
    if (~isempty(cst(n).logweighttable))
        % case of no parent: handled by never letting n == 1.
        parent = cst(cst(n).parent);

        % make sure ptable is initialized properly:
        % should have factors = cst factors
        if (isempty(parent.cpt.factors))
            error('initialize factors in discrete tables');
        end
        % load these values from the weighttable into the CondProbTable of
        % the parent:
        [vals, vind] = ValueIncrement(cst(n).dvalues, [], []);
        for i = 1:length(cst(n).logweighttable)
            % Stores weight in CondProbTable.factors so that we
            % know when one factor gets reduced to a single value by
            % evidence entry.  Otherwise the evidence gets lost in the CPT
            % and is not reflected in the ClusterSetTree:
            parent = parent.AddToCPT(cst(n).logweighttable(i), cst(n).dmembers, vind);

            [vals, vind] = ValueIncrement(cst(n).dvalues, vals, vind);
        end

        % select elements that correspond to the evidence provided by
        % EVVALS on the variables VARNAMES and factor-reduce this table
        % down to nodes present in the parent cluster:
        [parent.cpt, cst(n)] = factorreduceworker(parent.cpt, varnames, evvals, cst(n));

        % then take the prob. values from ptable and multiply them into the
        % weighttable of the parent
        
        % actually: don't need this info in the weighttable of the Cluster
        % Set Tree because it's in the Cond Prob Table of the CST.
        if (isempty(cst(cst(n).parent).logweighttable))
            cst(cst(n).parent).logweighttable = zeros(size(parent.cpt.logprob));
        end

        % save these back into the vars:
        cst(cst(n).parent) = parent;
    end    
end



