function cs = ForwardSimLearnParams(BN)
% vals = ForwardSim(BN)
% 
% generate a ClusterSetTree array for a BayesNet object for use with the
% ForwardSim() function.
%
% This is intended for use with Dynamic Bayes Nets that only have edges 
% from prior time points to future time points.
%
% INPUT :
%   BN : BayesNet object
%
% OUTPUT :
%   CS : ClusterSetTree array
%
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

    

% TO DO : 
% test with discrete/continuous mixed inference


%% for each node
nodes = BN.nodes;
cs = repmat(ClusterSetTree(),size(nodes));
for i = 1:length(nodes)
    n = nodes(i);
    
    %% make a ClusterSetTree for this node including the parents
    cs(i) = ClusterSetTree();
    cs(i).name = n.self;
    cs(i).index = n.index;
    neighbors = unique([n.pindex],'legacy');
    % a cluster is discrete if all members are discrete:
    cs(i).discrete = n.discrete;
    for j=1:length(neighbors)
        cs(i).members(end+1) = neighbors(j);
        if (nodes(neighbors(j)).discrete)
            % keep a separate list of discrete members:
            cs(i).dmembers(end+1) = neighbors(j);
            cs(i).dvalues{end+1} = nodes(neighbors(j)).values;
        end
    end

    % populate some CST values :
    cs(i) = InitClusterSetTree(cs(i), nodes);
    
    % use data just for this node and its parents. (Actually want all the
    % data for proper indexing)
    % initialize the CST using the data
    cs(i) = InitClusterPotential(cs(i),nodes,BN.data,BN.priorPrecision);
    cs(i) = PropagateCGRegression(cs(i), nodes);
    cs(i) = cs(i).InitStepSizes();
end
