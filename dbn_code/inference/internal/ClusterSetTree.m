classdef ClusterSetTree
% This class is a node in the junction tree.
% The junction TREE is an array of CLUSTERSET structs.
%
% MEMBERS : 
%         name % string
%         members % array of integers indexing into NODES of the bayes net
%         dmembers % array of discrete nodes in the cluster
%         index % index into NODES
%         parent % index into ClusterSetTree
%         children % array of indices into ClusterSetTree
%         discrete % boolean: true if ALL MEMBERS are discrete
%         dvalues % configurations of the discrete members that index the arrays:
%                 % a cell of arrays of actual discrete values each can
%                 % obtain
%         lppotential % cell of arrays of LPPOTs, one for each combo of discrete vars
%         postbag % another cell of arrays of LPPOTs
%         cpt % conditional prob table for discrete clusters
%         dvalstepsize % stepsizes for discrete value arrays
%                      % this represents a mixedbase integer, then take
%                      % dvalstepsize(i) is the PLACE value of the i^th digit.
%         activeflag % needed for PUSH evidence-entry algorithm
%         logweighttable % needed for storing potentials on boundary between
%                     % discrete and continuous nodes
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.



    properties
        name % string
        members % array of integers indexing into NODES of the bayes net
        dmembers % array of discrete nodes in the cluster
        index % index into NODES
        parent % index into ClusterSetTree
        children % array of indices into ClusterSetTree
        discrete % boolean: true if ALL MEMBERS are discrete
        dvalues % configurations of the discrete members that index the arrays:
                % a cell of arrays of actual discrete values each can
                % obtain
        lppotential % cell of arrays of LPPOTs, one for each combo of discrete vars
        postbag % another cell of arrays of LPPOTs
        cpt % conditional prob table for discrete clusters
        dvalstepsize % stepsizes for discrete value arrays
                     % this represents a mixedbase integer, then take
                     % dvalstepsize(i) is the PLACE value of the i^th digit.
        activeflag % needed for PUSH evidence-entry algorithm
        logweighttable % needed for storing potentials on boundary between
                    % discrete and continuous nodes
    end
    
    methods
        % constructor
        function cstree = ClusterSetTree(name, members, dmembers, index, parent, children, ...
                discrete, dvalues, lppotential, postbag, cpt, dsteps, activeflag, logweighttable)
            if (nargin > 0)
                cstree.name = name;
                cstree.members = members;
                cstree.dmembers = dmembers;
                cstree.index = index;
                cstree.parent = parent;
                cstree.children = children;
                cstree.discrete = discrete;
                cstree.dvalues = dvalues;
                cstree.lppotential = lppotential;
                cstree.postbag = postbag;
                cstree.cpt = cpt;
                cstree.dvalstepsize = dsteps;
                cstree.activeflag = activeflag;
                cstree.logweighttable = logweighttable;
            else
                cstree.name = '';
                cstree.members = [];
                cstree.dmembers = [];
                cstree.index = 0;
                cstree.parent = [];
                cstree.children = [];
                cstree.discrete = [];
                cstree.dvalues = [];
                cstree.lppotential = {[]};
                cstree.postbag = {[]};   
                cstree.cpt = CondProbTable();
                cstree.dvalstepsize  = [];
                cstree.activeflag = true;
                cstree.logweighttable = [];
            end
        end
            
        % include a method for converting the tree to a sparse adjacency
        % matrix
        function sparsemat = GetAdjacency(obj, treearray)
            % this function hardly needs to be part of this class/
            sparsemat = sparse(zeros(length(treearray)));
            for i = 1:length(treearray)
                for j = 1:length(treearray(i).children)
                    sparsemat(i,treearray(i).children(j)) = 1;
                end
                if (~isempty(treearray(i).parent))
                    sparsemat(treearray(i).parent, i) = 1;
                end
            end
        end
        
        
        % add an LPPotential object to each of the postbags of this CST:
        function obj = AddToLPPotential(obj, lppot)
            [obj,bag] = bagadder(obj, lppot, obj.lppotential);
            obj.lppotential = bag;
        end
        
        % add an LPPotential object to each of the postbags of this CST:
        function obj = AddToPostbags(obj, lppot)
            [obj,bag] = bagadder(obj, lppot, obj.postbag);
            obj.postbag = bag;
        end


        % add an LPPotential object to each of the postbags of this CST:
        function [obj, bagcell] = bagadder(obj, lppot, bagcell)
            % first see if there's a condition on this LPPot
            if (isempty(lppot.conditionvars))
                % if there's no condition, it goes in EVERY bag:
                if (~isempty(bagcell) && ~isempty(obj.dvalues))
                    [v, vinds] = ValueIncrement(obj.dvalues, []);
                    k = 1;
                    while (~isempty(v))
                        % lppot gets conditioned on the bag that it is put into:
                        rlppot = LPPotential(lppot.head, lppot.tail, lppot.params, ...
                            lppot.const, lppot.sigma, obj.dmembers, vinds);
                        [obj, bagcell] = bagadderworker(obj,rlppot, bagcell, k);
                        k = k+1;
                        [v, vinds] = ValueIncrement(obj.dvalues, v, vinds);
                    end
                else
                    [obj, bagcell] = bagadderworker(obj,lppot, bagcell, 1);
                end
            else
                % if there is a condition, find the bag(s) that have the
                % same condition:
                thiscondinds = zeros(size(lppot.conditionvars));
                thisvalinds = thiscondinds;
                noncondinds = true(size(obj.dmembers));
                for i = 1:length(lppot.conditionvars)
                    for j = 1:length(obj.dmembers)
                        if (obj.dmembers(j) == lppot.conditionvars(i))
                            thiscondinds(i) = j;
                            thisvalinds(i) = lppot.conditionvalinds(i);
                            noncondinds(j) = false;
                            break;
                        end
                    end
                end
                % now we have the indices of every input condition
                if (isempty(obj.dvalstepsize))
                    obj = obj.InitStepSizes();
                end
                increment = obj.dvalstepsize(thiscondinds) * (thisvalinds -1)';
                % loop through assignments to (dmembers) - (lppot.conditionvars)
                nonconds = (obj.dmembers(noncondinds));
                if (~isempty(nonconds))
                    bvec = zeros(size(nonconds));
                    numvalues = zeros(size(noncondinds));
                    for i = 1:length(noncondinds)
                        numvalues(i) = length(obj.dvalues{i});
                    end
                    basevec = numvalues(noncondinds);
                    while (~isempty(bvec))
                        % don't use basevec here, but the cumprod of it:
                        val = bvec * (obj.dvalstepsize(noncondinds))';
                        val = val + increment + 1; % val is computed for a zero-based array
                        % now, finally, add LPPot to the VAL^th postbag:
                        % need to EXPAND the CONDITION on these LPPots so
                        % to include everything that the current CST is
                        % conditioned upon:
                        rcond = [nonconds, lppot.conditionvars];
                        rcondvals = [bvec + 1, lppot.conditionvalinds];
                        rlppot = LPPotential(lppot.head, lppot.tail, lppot.params, lppot.const, ...
                            lppot.sigma, rcond, rcondvals);
                        [obj, bagcell] = bagadderworker(obj,rlppot, bagcell, val);
                        bvec = bitincbasevec(bvec,basevec);
                    end
                else
                    % all conditioning variables are present, so the
                    % condition in the LPPot defines a unique bagcell
                    ind = increment + 1;
                    [obj, bagcell] = bagadderworker(obj,lppot, bagcell, ind);
                end
            end
        end
        
        function [obj, bagcell] = bagadderworker(obj, lppot, bagcell, ind)
            if (isempty(bagcell))
                bagcell = {[lppot]};
            elseif (isempty(bagcell{ind}))
                bagcell{ind} = [lppot];
            else
                bagcell{ind}(end+1) = lppot;
            end 
        end
        
        
        % similar to the routine to add an LPPot to the bags, but do it
        % with adding a logweighttable to the CondProbTable object of the CST.
        function obj = AddToCPT(obj, logweight, condvars, condvals)
            % input: 
            % condvars are the calling (different from obj) CST.dmembers
            % condvals are indices into calling CST.dvalues.  such that
            % condvars(i) takes values cst.dvalues{find(cst.dmembers ==
            % condvars(i))}(condvals(i)).
            % first see if there's a condition on this LPPot
            if (isempty(condvars))
                % if there's no condition, it goes in EVERY bag:
                if (~isempty(obj.cpt.logprob))
                    [v, vinds] = ValueIncrement(obj.cpt.values, []);
                    k = 1;
                    while (~isempty(v))
                        % find the right entry in the cpt.logprob:
                        obj.cpt.logprob(k) = obj.cpt.logprob(k) + logweight;
                        k = k+1;
                        [v, vinds] = ValueIncrement(obj.cpt.values, v, vinds);
                    end
                else
                    error('CPT had no probability table');
                end
            else
                % if there is a condition, match those conditions
                thiscondinds = zeros(size(condvars));
                thisvalinds = thiscondinds;
                noncondinds = true(size(obj.cpt.findex));
                for i = 1:length(condvars)
                    for j = 1:length(obj.cpt.findex)
                        if (obj.cpt.findex(j) == condvars(i))
                            thiscondinds(i) = j;
                            thisvalinds(i) = condvars(i);
                            noncondinds(j) = false;
                            break;
                        end
                    end
                end
                allindices = 1:length(obj.cpt.findex);
                noncondinds = allindices(noncondinds);
                % now we have the indices of every input condition
                obj.cpt = obj.cpt.InitStepSizes();
                vallengths = ones(size(obj.cpt.values));
                for i = 1:length(obj.cpt.values)
                    vallengths(i) = length(obj.cpt.values{i});
                end
                [allinds] = allvaluesindex(noncondinds, vallengths, obj.cpt.stepsize,...
                    thiscondinds, condvals);
                % maybe thiscondinds(thiscondinds > 0) ? 
                obj.cpt.logprob(allinds) = obj.cpt.logprob(allinds) + logweight;
            end
            
        end
 
        % enter evidence that DISCVARS have been assigned values DISCVALS
        function obj = DiscConditionEvidence(obj, discvars, discvals)
            % just narrow down everything to only have the remaining
            % conditions
            
            % just return if the discvars don't overlap with the discrete
            % conditioning vars in this factor:
            if (isempty(intersect(obj.dmembers, discvars)))
                return;
            end

            % now we have the indices of every input condition
            if (isempty(obj.dvalstepsize))
                obj = obj.InitStepSizes();
            end
            
            % find the local index of discvars
            thiscondinds = zeros(size(discvars));
            thisvalinds = thiscondinds;
            noncondinds = true(size(obj.dmembers));
            for i = 1:length(discvars)
                for j = 1:length(obj.dmembers)
                    if (obj.dmembers(j) == discvars(i))
                        thiscondinds(i) = j;
                        for k = 1:length(obj.dvalues{j})
                            if (discvals(i) == obj.dvalues{j}{k})
                                thisvalinds(i) = k;
                                break;
                            end
                        end
                        noncondinds(j) = false;
                        break;
                    end
                end
            end
            
            % find the discvalue index of discvals
            
            % discvars unassigned on this factor
            allinds = 1:length(obj.dmembers);
            unassingedvars = allinds(noncondinds);
            % this gives all the indices into mixedbase linearized arrays
            % that have the unassigned vars unassinged and the assigned
            % vars given values by thisvalinds
            vallengths = ones(size(obj.dvalues));
            for i = 1:length(obj.dvalues)
                vallengths(i) = length(obj.dvalues{i});
            end
            % sort these incase they're mis-ordered:
            thiscondinds = thiscondinds(thiscondinds > 0);
            thisvalinds = thisvalinds(thisvalinds > 0);
            [thiscondinds,smap] = sort(thiscondinds);
            thisvalinds = thisvalinds(smap);
            keepinds = allvaluesindex(unassingedvars, vallengths, obj.dvalstepsize,...
                thiscondinds, thisvalinds);
            
            % now just drop everything that's in the dropinds:
            if (max(keepinds) > length(obj.lppotential))
                error('Discrete Evidence conditioning in Continuous CST is greater than existing LPPotentials');
            end
            obj.lppotential = obj.lppotential(keepinds);
            obj.postbag = obj.postbag(keepinds);
            if (~isempty(obj.logweighttable))
                obj.logweighttable = obj.logweighttable(keepinds);
            end        
        end
        
        %compute stepsize for each of these
        function obj = InitStepSizes(obj)
            obj.dvalstepsize(1) = 1;
            for i = 2:length(obj.dmembers)
                obj.dvalstepsize(i) = obj.dvalstepsize(i-1) * length(obj.dvalues{i-1}); 
            end
        end
            
    end
        
end
