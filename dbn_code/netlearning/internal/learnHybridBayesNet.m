function [BN,outstats]=learnHybridBayesNet(contData, discData, priorPrecision, ...
         contphencol, searchParameter, initialBayesNet, verbose)
% [BN,outstats]=learnHybridBayesNet(contData,discData,priorPrecision,searchParameter)
% learns the optimal Bayesian network given a continuous data set and a 
% discrete data set. The Bayesian network learning is implemented by 
% stepwise K2-algorithm.
%
% INPUT:
% CONTDATA: numeric data matrix of continous nodes (indexed by node,sample)
% DISCDATA: numeric data matrix of discrete nodes (indexed by node,sample)
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% CONTPHENCOL: if we have a continuous phenotype, this must be the index of
%   that (row) into the the CONTDATA array.  Otherwise, set CONTPHENCOL =
%   [];
% SEARCHPARAMETER:
%   searchParameter.stepwise: true=use stepwise search; false=not use
%   searchParameter.MAX_NUM_PARENTS: maximum number of parents per node
%   searchParameter.BFTHRESH: minimum increase in Log Likelihood for 
%       inclusion of an edge in the network.  Deafult = 0;
% INITIALBAYESNET: initial network structure, to continue the search for
%   additonal edges on an existing network.
% VERBOSE: if true, increases output.
%
% OUTPUT:
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

if (nargin < 7)
    verbose = false;
end

%% check input arguments
if nargin < 6
    error('Syntax error!');
    return;
end

if (~isfield(searchParameter,'BFTHRESH'))
    searchParameter.BFTHRESH = 0;
end


%% check if both types of data are present
if isempty(contData) && isempty(discData)
    error('No data available for learning hybrid Bayesian network!');
    return;
end




%% data info
numContNode = size(contData,1);
numDiscNode = size(discData,1);
numAllNode = numContNode+numDiscNode;
numSampleC = size(contData,2);
numSampleD = size(discData,2);
if ~isempty(contData) && ~isempty(discData) && numSampleC~=numSampleD
    error('Sample sizes of the continuous and discrete data must be equal!');
end
clear numSample*;






%% preallocate for storing results
% note: adjMatrix(m,n)=1 if m is a parent of n; 
% note: adjMatrix: the columns/rows are indexed by continuous nodes followed by discrete nodes
%
if ~isempty(initialBayesNet) %if previous network structure is given
    adjMatrix = initialBayesNet.adjMatrix;
    nodeModel = initialBayesNet.nodeModel;
    contParent= initialBayesNet.contParent;
    discParent= initialBayesNet.discParent;
    bothParent= initialBayesNet.bothParent;        
    edgeweight= initialBayesNet.edgeweight;        
    weightMatrix= initialBayesNet.weightMatrix;        
else    
    adjMatrix = spalloc(numAllNode,numAllNode,numAllNode*searchParameter.MAX_NUM_PARENTS);
    weightMatrix= spalloc(numAllNode,numAllNode,numAllNode*searchParameter.MAX_NUM_PARENTS);
    nodeModel = cell(1,numAllNode);
    contParent= cell(1,numAllNode);
    discParent= cell(1,numAllNode);
    bothParent= cell(1,numAllNode);
    edgeweight= cell(1,numAllNode);
end

% search diagnostic parameters
lldiffs = zeros(1,numContNode+1);
evals = zeros(1,numContNode+1);
nedges = zeros(1,numContNode+1);
addededges = cell(1,numContNode+1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start learning network of continuous nodes 
% (note: continuous nodes are allowed to have discrete and continuous parents)

%% search best parent(s) for each continuous node 
for g=1:numContNode  %note: this loop is skipped if no continuous data 
    if (verbose)
        fprintf(1,'\t learning parent(s) of continuous node %d\n',g);
    end
    
    %% data of the child node
    childData = contData(g,:);       

    
    %% initial model of node g 
    if isempty(contParent{g})
        contParentData=[];
    else
        contParentData=contData(contParent{g},:);        
    end
    if isempty(discParent{g})
        discParentData=[];
    else
        discParentData=discData(discParent{g},:);        
    end
    model_ini=learnLocalBN_MixToCont(contParentData,discParentData,childData,priorPrecision);
    
    %% search 
    if (~isempty(contphencol) && g == contphencol)
        done = true;
    else
        done = false;
    end
    numevals = 0;
    lltempdiff = 0;
    while ~done && length(bothParent{g})<searchParameter.MAX_NUM_PARENTS
        allPossibleParent = g+1:numAllNode;
        %candidate nodes: potential parents to be added in BN model
        candidateNode=setdiff(allPossibleParent,bothParent{g},'legacy');

        if ~isempty(candidateNode)

            %% add each candidate as a parent of the current child node
            model_addnew = cell(1,length(candidateNode));
            logLLH_addnew=zeros(1,length(candidateNode));
            for a=1:length(candidateNode)
                bothParent_addnew = [bothParent{g},candidateNode(a)];                
                contParentData_addnew=contData(bothParent_addnew(bothParent_addnew<=numContNode),:);
                discParentData_addnew=discData(bothParent_addnew(bothParent_addnew>numContNode)-numContNode,:);
                model_addnew{a}=learnLocalBN_MixToCont(...
                                contParentData_addnew,discParentData_addnew,childData,priorPrecision);
                logLLH_addnew(a)=model_addnew{a}.logLLH;
            end
            numevals = numevals + length(candidateNode);

            %% choose which added parent resulting max log-likelihood
            [maxLLH, maxInd_addnew]=max(logLLH_addnew);
            maxInd_addnew=maxInd_addnew(1); %in case of multiple maxima
            newAddParent=candidateNode(maxInd_addnew); 
            
            
            %% tentatively treat adding newAddParent as the best model
            bothParent_best=[bothParent{g}, newAddParent];
            % also add edge weight here:
            edgeweight_best=[edgeweight{g}, maxLLH - model_ini.logLLH];
            model_best = model_addnew{maxInd_addnew};
            clear *_addnew;
            
            
            
            %% stepwise search
            if searchParameter.STEPWISE && (length(bothParent{g})>1) 
                
                %% preallocate
                bothParent_rem=cell(1,length(bothParent{g}));
                logLLH_rem=zeros(1,length(bothParent{g}));
                model_rem=cell(1,length(bothParent{g}));
                edgeweight_rem=cell(1,length(bothParent{g}));
                
                %% remove a parent, add a potential parent, calculate logLLH
                for a=1:length(bothParent{g})
                    % remove a parent, add the potential parent
                    bothParent_rem{a}=bothParent_best;
                    bothParent_rem{a}(a)=[];
                    edgeweight_rem{a}=edgeweight_best;
                    edgeweight_rem{a}(a)=[];

                    % parent data
                    contParentData_rem=contData(bothParent_rem{a}(bothParent_rem{a}<=numContNode),:);
                    discParentData_rem=discData(bothParent_rem{a}(bothParent_rem{a}>numContNode)-numContNode,:);
                    model_rem{a}=learnLocalBN_MixToCont(...
                                  contParentData_rem,discParentData_rem,childData,priorPrecision);                    
                    logLLH_rem(a)=model_rem{a}.logLLH;
                end
                numevals = numevals + length(bothParent{g});

                %% check if no-parent-removal is better
                [maxLogLLH_rem, maxInd_rem]=max(logLLH_rem);
                maxLogLLH_rem=maxLogLLH_rem(1);
                maxInd_rem=maxInd_rem(1);
                if (maxLogLLH_rem>=model_best.logLLH)
                    bothParent_best=bothParent_rem{maxInd_rem};
                    edgeweight_best=edgeweight_rem{maxInd_rem};
                    edgeweight_best(end) = maxLogLLH_rem;
                    model_best = model_rem{maxInd_rem};
                end 
                clear *_rem;
                
            end %end stepwise search

            
            %% check if model_best better than the initial model 
            if (model_best.logLLH - searchParameter.BFTHRESH) > model_ini.logLLH
                % with this edge
                lltempdiff = lltempdiff + model_best.logLLH - model_ini.logLLH;
                model_ini = model_best;
                bothParent{g} = bothParent_best;
                contParent{g} = bothParent_best(bothParent_best<=numContNode);
                discParent{g} = bothParent_best(bothParent_best>numContNode)-numContNode;
                edgeweight{g} = edgeweight_best;
            else
                done=1;
            end
            clear *_best;

        else
            done=1;  %if there are no more potential parents, exit cycle
        end %end IF ~isempty(candidateNode)

    end %end WHILE (cycle on node parents)


    % update results
    nodeModel{g} = model_ini;
    if (~isempty(bothParent{g}))
        % bothParent{g} can be an array, but it doesn't matter here since
        % a scalar assignment to a matrix works fine
        adjMatrix(bothParent{g},g)=1;
        % this assingment is an array to an array:
        weightMatrix(bothParent{g},g) = edgeweight{g};
    end
    nedges(g+1) = nedges(g) + length(bothParent{g});
    evals(g+1) = evals(g) + numevals;
    lldiffs(g+1) = lldiffs(g) + lltempdiff;
    edges = cell(1,length(bothParent{g}));
    for c = 1:length(bothParent{g})
        edges{c} = [bothParent{g}(c),g];
    end
    addededges{g+1} = edges;
    clear model_ini;
    
end %end FOR (loop on nodes)

%% end learning network of continuous nodes 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% learn network of discrete nodes (note: discrete nodes must not have continuous 
%% parents, so we just need to use the function that handles discrete nodes)
if numDiscNode>0
    [BayesNetDisc, discstats]=learnDiscreteBayesNet(discData, priorPrecision, searchParameter);

    %% combine networks of continuous and discrete nodes
    adjMatrix(numContNode+1:end,numContNode+1:end)=BayesNetDisc.adjMatrix;
    weightMatrix(numContNode+1:end,numContNode+1:end)=BayesNetDisc.weightMatrix;
    nodeModel(numContNode+1:end) = BayesNetDisc.nodeModel;
    discParent(numContNode+1:end) = BayesNetDisc.discParent;
    bothParent(numContNode+1:end) = discParent(numContNode+1:end);
    
    lldiffs = [lldiffs, discstats.lldiffs + lldiffs(end)];
    evals = [evals, discstats.numevals + evals(end)];
    nedges = [nedges, discstats.numedges + nedges(end)];
    addededges = {addededges{:}, discstats.addededges{:}};
end


%% output
BN.adjMatrix = adjMatrix;
BN.nodeModel = nodeModel;
BN.discParent = discParent;
BN.contParent = contParent;
BN.bothParent = bothParent;
BN.weightMatrix = weightMatrix;

outstats.lldiffs = lldiffs;
outstats.numevals = evals;
outstats.numedges = nedges;
outstats.addededge = addededges;
return;