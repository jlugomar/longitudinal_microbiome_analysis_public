classdef CondProbTable
% Defines the CONDPROBTABLE class
%
% MEMBERS:
% cpt.node : cell string of main node name
% cpt.factors : cell array of strings of all node names included
%               factors{1} = node{1}
% cpt.table : 2-dim array of values and factor indices
%             table(:,i) are values of factor{i}  
% cpt.logprob : array of conditional log probabilities such that if all factors
%            take on values table(k,:) then factors{1} has log probability 
%            logprob(k).  this is unnormalized.
% cpt.values : cell array of arrays, where values{i} is an array of
%           values that factors{i} assumes
% cpt.findex : indices into NODES for each FACTOR included
% cpt.stepsize : internal variable for keeping track of indexing, equal to
%           the cumprod of lengths(cpt.values).  Array of numbers.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

    properties
        node 
        factors
        table
        logprob
        values
        findex
        stepsize
    end
    
    methods
        % constructor
        function cpt = CondProbTable(node, factors, table, logprob, values, ...
                findex, stepsize)
            if (nargin == 0)
                cpt.node = '';
                cpt.factors = {};
                cpt.table = [];
                cpt.logprob = [];
                cpt.values = {};
                cpt.findex = [];
                cpt.stepsize = [];
            else
                cpt.node = node;
                cpt.factors = factors;
                cpt.table = table;
                cpt.logprob = logprob;
                cpt.values = values;
                cpt.findex = findex;
                cpt.stepsize = stepsize;
            end
        end
        
        %compute stepsize for each of these
        function obj = InitStepSizes(obj)
            obj.stepsize(1) = 1;
            for i = 2:length(obj.findex)
                obj.stepsize(i) = obj.stepsize(i-1) * length(obj.values{i-1}); 
            end
        end

    end
    
end
