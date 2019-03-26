function cstnode = InitContFactor(cstnode, node, nodes, data, priorPrecision)
%cstnode = InitContFactor(cstnode, node, nodes, data, priorPrecision)
%
% Initializes a continuous or mixed discrete/continuous factor for a
% ClusterSetTree item.
% 
% INPUT : 
%   CSTNODE : an instantiation of the ClusterSetTree class, this represents
%      the node in the junction tree to be initialized
%   NODE : the node corresponding to the CSTNODE
%   NODES : the master list of nodes in the bayesian network
%   DATA : an array of values for each variable in the domain of the bayes
%       net
%   PRIORPRECISION : hyperparameters for Bayesian priors
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% find each member of the factor in the datafile - a column
memcols = zeros(size(node.parents));
discols = [];  % these are indices of memcols, which is reduced from DATA
contcols = [];
contmembers = [];
disvalues = [];
dismembers = [];
for i = 1:length(node.parents)
    memcols(i) = nodes(node.pindex(i)).colind;
    if (nodes(node.pindex(i)).discrete)
        discols(end+1) = i;
        disvalues{end+1} = nodes(node.pindex(i)).values;
        dismembers(end+1) = node.pindex(i);
    else
        contcols(end+1) = i;
        % list of continuous dependecies:
        contmembers(end+1) = node.pindex(i);
    end
end

% "members" does not contain the self index
selfdata = data(:,node.colind);
data = data(:, memcols);

[v, vinds, numvalcombos] = ValueIncrement(disvalues, []);
% update hyperparemeter "nu" in priorPrecision
priorPrecision.nu = priorPrecision.nu/numvalcombos;
for d = 1:numvalcombos
    discdata = data(:,discols);
    vmat = repmat(v, size(discdata,1),1);
    if (isempty(vmat))
        % this is the case of zero discrete parents:
        zs = zeros(size(selfdata));
    else
        matches  = discdata - vmat;
        % matches are those rows who's dot product is zero with itself
        zs = diag(abs(matches) * abs(matches'));
    end
    % now we can take the matching rows and do a REGRESSION on them
    if (sum(zs == 0) > 0)
        y = selfdata(zs == 0);
        % contmembers is a subset of the PARENTS of the NODE.
        if (~isempty(contmembers))
            xs = data(zs == 0, contcols);
            % call bayesian prior function; this requires data to be
            % transposed:            
            localModel = learnLocalBN_ContToCont(xs',y',priorPrecision);
            b = localModel.regreCoeff;
            if (length(b) < 2)
               error('suspected dependence unfound!  probably dependence is close to zero!');
            end

            % use a biased estimate of variance to handle the case of low data
            %var_y = localModel.postVar;
            var_y = 1/localModel.postPrec;
            if (var_y < 0)
                error('Error: negative variance estimate from prior!');
            end
            
            xpot = LPPotential(node.index, contmembers, b(2:end), b(1), var_y, ...
                dismembers, vinds);
        else
            % still use bayesian priors:
            localModel = learnLocalBN_ContToCont([],y',priorPrecision);
            b = localModel.regreCoeff;
            % use a biased estimate of variance to handle the case of low data
            %var_y = localModel.postVar;
            var_y = 1/localModel.postPrec;
            if (var_y < 0)
                error('Error: negative variance estimate from prior!');
            end
            % ok to extend the conditioning to cstnode.dmembers here
            xpot = LPPotential(node.index, [], [], b, var_y, ...
                dismembers, vinds);
        end
    else
        % no data fitting the configuration of the discrete nodes -
        % uses the prior from priorPrecision to compute the base probability 
        % according to the prior, w/o any data:
        
        % this just uses the global mean as an estimate of the y-value here
        localModel = learnLocalBN_ContToCont([],mean(selfdata),priorPrecision);
        b = localModel.regreCoeff;
        %var_y = localModel.postVar;
        var_y = 1/localModel.postPrec;
        if (var_y < 0)
            error('Error: negative variance estimate from prior!');
        end
        if (length(b)-1 < length(contmembers))
            b = [b;zeros(length(contmembers)- (length(b)-1),1)];
        end
        xpot = LPPotential(node.index, contmembers, b(2:end), b(1), var_y, ...
            dismembers, vinds);
    end
    % now check to see if this cluster's ELIMINATION NODE is x:
    if (cstnode.index == node.index)
        cstnode = cstnode.AddToLPPotential(xpot);
    else
        cstnode = cstnode.AddToPostbags(xpot);
    end
    [v, vinds] = ValueIncrement(disvalues, v, vinds);
end
