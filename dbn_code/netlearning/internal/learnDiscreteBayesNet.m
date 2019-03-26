function [BN, outstats]=learnDiscreteBayesNet(discData, priorPrecision, searchParameter, verbose)
% [BN, outstats]=learnDiscreteBayesNet(discData, priorPrecision, searchParameter, verbose)
% learns the optimal Bayesian network given a discrete data set.
% The Bayesian network learning is implemented by K2-algorithm.
%
% INTERNAL
%
% Input:
% DISCDATA: numeric data matrix of discrete nodes (indexed by node,sample)
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% SEARCHPARAMETER:
%   searchParameter.stepwise: true=use stepwise search; false=not use
%   searchParameter.MAX_NUM_PARENTS: maximum number of parents per node
%   searchParameter.BFTHRESH: minimum increase in Log Likelihood for 
%       inclusion of an edge in the network.  Deafult = 0;
% VERBOSE: if true, increases output.
%
% Output:
% BN: not a BayesNet class object; but rather a structure with fields: 
%	BN.adjMatrix: network structure recorded by adjacency matrix:
%   	adjMatrix(m,n)=1 if m is a parent of n.
%   BN.nodeModel: array of model configuration structures, one per
%   	node. Includes field BN.nodeModel.logLLH
%	BN.weightMatrix: network structure recorded by weighted adjacency matrix:
%   	weightMatrix(m,n)=X, where X is the change in log likelihood
%   	associated with adding that edge to the network; considering the other
%       edges that were present at the time that edge was added to the network.
% OUTSTATS: structure with fields describing the characteristics of the
%       search procedure, in arrays per "step;" a step is either an edge
%       added or removed.
%   outstats.lldiffs: difference in loglikelihood at each step of the algorithm.
%   outstats.numedges: number of edges in the network at each step of the algorithm.
%   outstats.numevals: number of potential network states evaluated at each
%       step.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.


if (nargin < 4)
    verbose = false;
end

if (~isfield(searchParameter,'BFTHRESH'))
    searchParameter.BFTHRESH = 0;
end


%% check if data is present
if isempty(discData) 
    error('No data available for learning discrete Bayesian network!');
    outstats.lldiffs = [0];
    outstats.numevals = [0];
    outstats.numedges = [0];
    outstats.addededges = {};
    return;
end



%% check if more than one node in the data
if  size(discData,1)==1
    BN.adjMatrix(1,1) = 0;
    BN.nodeModel{1} = [];    
    BN.discParent{1} = [];
    BN.weightMatrix(1,1) = 0;
    outstats.lldiffs = [0];
    outstats.numevals = [0];
    outstats.numedges = [0];
    outstats.addededges = {};
    return;
end



%% data info
numDiscNode = size(discData,1);



%% preallocate for storing results
adjMatrix = spalloc(numDiscNode,numDiscNode,numDiscNode*searchParameter.MAX_NUM_PARENTS);
weightMatrix = spalloc(numDiscNode,numDiscNode,numDiscNode*searchParameter.MAX_NUM_PARENTS);
nodeModel = cell(1,numDiscNode);
discParent=cell(1,numDiscNode);
edgeweight= cell(1,numDiscNode);

% search diagnostic parameters
lldiffs = zeros(1,numDiscNode);
evals = zeros(1,numDiscNode);
nedges = zeros(1,numDiscNode);
addededges = cell(1,numDiscNode);

%% search best parent(s) for each node 
for g=1:numDiscNode-1
    if (verbose)
        fprintf(1,'\t learning parent(s) of discrete node %d\n',g);
    end
    
    %% data of the child node
    childData = discData(g,:);        

    
    %% initial model of node g (i.e., no parent)
    discParent{g}=[];
    discParentData=[];    
    model_ini=learnLocalBN_DiscToDisc(discParentData,childData,priorPrecision);
    
        
    
    %% search 
    done = 0; 
    numevals = 0;
    lltempdiff = 0;
    while ~done && length(discParent{g})<searchParameter.MAX_NUM_PARENTS
        
        allPossibleParent = g+1:numDiscNode;
        %candidateDiscNode: potential discrete parents to be added in BN model
        candidateDiscNode=setdiff(allPossibleParent,discParent{g},'legacy');
        
        numevals = 0;
        if ~isempty(candidateDiscNode)

            %% add each candidate as a parent of the current child node
            model_addnew = cell(1,length(candidateDiscNode));
            logLLH_addnew=zeros(1,length(candidateDiscNode));
            for a=1:length(candidateDiscNode)
                discParent_addnew = [discParent{g},candidateDiscNode(a)];
                model_addnew{a}=learnLocalBN_DiscToDisc(discData(discParent_addnew,:),childData,priorPrecision);
                logLLH_addnew(a)=model_addnew{a}.logLLH;
            end
            numevals = numevals + length(candidateDiscNode);

            %% choose which added parent resulting max log-likelihood
            [maxLLH, maxInd_addnew]=max(logLLH_addnew);
            maxInd_addnew=maxInd_addnew(1); %in case of multiple maxima
            newAddParent=candidateDiscNode(maxInd_addnew); 
            
            
            %% tentatively treat adding newAddParent as the best model
            discParent_best=[discParent{g}, newAddParent];
            model_best = model_addnew{maxInd_addnew};
            % also add edge weight here:
            edgeweight_best=[edgeweight{g}, maxLLH - model_ini.logLLH];
            clear *_addnew;
            
            
            
            %% stepwise search
            if searchParameter.STEPWISE && (length(discParent{g})>1) 
                
                %% preallocate
                discParent_rem=cell(1,length(discParent{g}));
                discLogLLH_rem=zeros(1,length(discParent{g}));
                discModel_rem=cell(1,length(discParent{g}));
                edgeweight_rem=cell(1,length(discParent{g}));
                
                %% remove a parent, add a potential parent, calculate logLLH
                for a=1:length(discParent{g})
                    % remove a parent, add the potential parent
                    discParent_rem{a}=discParent_best;
                    discParent_rem{a}(a)=[];
                    edgeweight_rem{a}=edgeweight_best;
                    edgeweight_rem{a}(a)=[];

                    % parent data: remove the data of parent a
                    discModel_rem{a}=learnLocalBN_DiscToDisc(discData(discParent_rem{a},:),childData,priorPrecision);
                    discLogLLH_rem(a)=discModel_rem{a}.logLLH;
                end
                numevals = numevals + length(discParent{g});

                %% check if no-parent-removal is better
                [maxLogLLH_rem, maxInd_rem]=max(discLogLLH_rem);
                maxLogLLH_rem=maxLogLLH_rem(1);
                maxInd_rem=maxInd_rem(1);
                if (maxLogLLH_rem>=model_best.logLLH)
                    discParent_best=discParent_rem{maxInd_rem};
                    edgeweight_best=edgeweight_rem{maxInd_rem};
                    edgeweight_best(end) = maxLogLLH_rem;
                    model_best = discModel_rem{maxInd_rem};
                end 
                clear *_rem;
                
            end %end stepwise search

            
            %% check if model_best better than the initial model 
            if (model_best.logLLH - searchParameter.BFTHRESH) >model_ini.logLLH
                lltempdiff = lltempdiff + model_best.logLLH - model_ini.logLLH;
                discParent{g} = discParent_best;
                edgeweight{g} = edgeweight_best;
                model_ini = model_best;
            else
                done=1;
            end
            clear *_best;

        else
            done=1;  %if there are no more potential parents, exit cycle
        end %end IF ~isempty(candidateDiscNode)

    end %end WHILE (cycle on node parents)


    % update results
    nodeModel{g} = model_ini;
    if (~isempty(discParent{g}))        
        adjMatrix(discParent{g},g)=1;
        weightMatrix(discParent{g},g) = edgeweight{g};
    end
    edges = cell(1,length(discParent{g}));
    for c = 1:length(discParent{g})
        edges{c} = [discParent{g}(c),g];
    end
    addededges{g+1} = edges;
    nedges(g+1) = nedges(g) + length(discParent{g});
    evals(g+1) = evals(g) + numevals;
    lldiffs(g+1) = lldiffs(g) + lltempdiff;
end %end FOR (loop on nodes)



%% output
BN.adjMatrix = adjMatrix;
BN.nodeModel = nodeModel;
BN.discParent = discParent;
BN.weightMatrix = weightMatrix;

outstats.lldiffs = lldiffs;
outstats.numevals = evals;
outstats.numedges = nedges;
outstats.addededges = addededges;

return;