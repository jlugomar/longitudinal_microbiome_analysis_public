function [BayesNet, outstats] = RandomizeHillClimbingNetSearchGlobal(contData, discData, ...
    priorPrecision, phencol, searchParameter)
%[BayesNet, outstats] = ExhaustiveFullNetSearch(contData, discData, priorPrecision, phencol, searchParameter)
%
% Performs a hillclimbing search over all possible edges, adding one 
% edges at a time.  In contrast to the usual learning algorithm
% that uses a K2-style node order, this one builds the most likely
% network across all input variables.  Allows backtracking, when removing
% an edge results in an increase in posterior log likelihood of the data.  
% Roughly O(k(n^2)*d*f), where k = priorPrecision.maxParents, n = number
% of nodes, d = size of datasample, and f is the time required for DFS for
% cycle checking, approx ~ O(log(n)).
%
% Input:
% CONTDATA: numeric data matrix of continous nodes (indexed by node,sample)
% DISCDATA: numeric data matrix of discrete nodes (indexed by node,sample)
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents for each node
%       priorPrecision.mle; % if true, maximum likelihood estimation will be used
% PHENCOL: index indicating the phenotype column.  Either the last contdata
%   column or the last discdata column (with contdata before discdata)
% SEARCHPARAMETER: structure with the following fields:
%   searchParameter.BF_THRESH: a Bayes Factor threshold for including an 
%       edge.  Edges with likelihood below BF_THRESH will not be included.  
%       This is the main stopping criterion.
%   searchParameter.nophenotype: a boolean, if true, means that the
%       network should be built without any particular phenotype
%   searchParameter.backtracking: a boolean, if true, means that the
%       network should be built allowing backtracking.  Default = true.
%   searchParameter.DBN: a boolean, if true, indicates that we're learing a
%       two-stage dynamic bayes net, and will use the TSDBNStateTrackerSearch
%       class. Default = false.
%       A dynamic bayes net must further define:
%   searchParameter.t0cont: continuous vars of the first of two-stages
%   searchParameter.tncont: continuous vars of the second of two-stages
%   searchParameter.t0disc: discrete vars of the first of two-stages
%   searchParameter.tndisc: discrete vars of the second of two-stages        
%
% OUTPUT:
% BAYESNET: structure with field:
%   BayesNet.adjMatrix: adjacency matrix for the bayes net.    
% OUTSTATS: structure with fields describing the characteristics of the
%       search procedure, in arrays per "step;" a step is either an edge
%       added or removed.
%   outstats.lldiffs: difference in loglikelihood at each step of the algorithm.
%   outstats.numedges: number of edges in the network at each step of the algorithm.
%   outstats.numevals: number of potential network states evaluated at each
%       step.
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

% Phenotype node cannot have any edges going in
if (~isfield(searchParameter,'nophenotype'))
    nophenotype = false;
else
    nophenotype = searchParameter.nophenotype;
end

if (~isfield(searchParameter,'backtracking'))
    backtracking = true;
else
    backtracking = searchParameter.backtracking;
end

if (~isfield(searchParameter,'BF_THRESH'))
    bf_thresh = 0;
    searchParameter.BF_THRESH = 0;
else
    bf_thresh = searchParameter.BF_THRESH;
end

if (~isfield(searchParameter,'DBN'))
    DBN = false;
else
    DBN = searchParameter.DBN;
end

if (~isfield(searchParameter,'unwrapped'))
    unwrapped = false;
else
    unwrapped = searchParameter.unwrapped;
end

% Set up the network structure search
search = StateTrackerSearch(contData, discData, priorPrecision, phencol, backtracking, bf_thresh, nophenotype);
if (DBN)
    if (unwrapped)
        search = UnwrappedTSDBNSTS(contData, discData, priorPrecision, phencol, backtracking, bf_thresh, nophenotype);
        search = search.Init();
        search = search.SetUnwrappedTimeSeries();
    else
        % A Dynamic BN must define:
        %   searchParameter.t0cont
        %   searchParameter.tncont
        %   searchParameter.t0disc
        %   searchParameter.tndisc
        
        % Set up the search to do a dynamic bayesian network
        search = TSDBNStateTrackerSearch(contData, discData, priorPrecision, phencol, backtracking, bf_thresh, nophenotype);
        search = search.Init();
        search = search.SetTimeSeries(searchParameter.t0cont, searchParameter.t0disc, ...
            searchParameter.tncont, searchParameter.tndisc);
    end
else
    search = search.Init();
end

verbose = false;
llscores = 0;
numedges = 0;
addedeges = {};
done = false;
invalidOperation = false;
evalcount = 0;
evalnums = 0;
maxScore = -Inf;
currentScore = 0;
numSample = 0.5 * (size(contData, 2) + size(discData, 2));

%% Initialize network structure with self-loops (only for DBNs)
if (isfield(searchParameter, 'DBN') && searchParameter.DBN)
    % Check all possible self-edges
    [search, numevals] = search.EvalSelfEdges();
    evalcount = evalcount + numevals;
    candidateSelfEdges = find(search.lldiff > bf_thresh);
    % If no more candidate self-edges, exit initialization
    if (isempty(candidateSelfEdges))
        done = true;
    end
    % Check each of these candidate self-edges
    i = 0;
    while (~done && i < length(candidateSelfEdges))
        i = i + 1;
        currentSelfEdge = candidateSelfEdges(i);
        [choseni, chosenj] = ind2sub(size(search.llscore), currentSelfEdge);
        currentEdgeScore = search.llscore(choseni, chosenj);
        currentEdgeDiffScore = search.lldiff(choseni, chosenj);
        oldEdgeScore = currentEdgeScore - currentEdgeDiffScore;
        % Decide if a candidate edge should be added. If so, add edge.
        [search, success] = search.DoEdgeAdd(choseni, chosenj, currentEdgeDiffScore);
        if (success)
            % Debugging
            if (verbose)
                outline = ['Edge op succesful;', num2str(choseni), ';', num2str(chosenj), ';', num2str(currentEdgeDiffScore)];
                disp(outline);
            end
            currentScore = llscores(end) + currentEdgeScore;
            llscores(end + 1) = llscores(end) + currentEdgeScore;
            addedeges{end + 1} = [choseni, chosenj];
            numedges(end + 1) = search.edges;
            evalnums(end + 1) = evalcount;
        end
    end
end

% Compute number of free parameters for candidate network structure
[search, currentNumParameters] = search.GetNumParameters(choseni, chosenj, false);
% Compute BIC penalty score
currentPenaltyScore = (0.5 * currentNumParameters) * log(numSample);
if (~searchParameter.bic)
    % Compute AIC penalty score
    currentPenaltyScore = currentNumParameters;
end
% Compute initial network structure BIC or AIC score 
currentPenalizedScore = currentScore - currentPenaltyScore;

%% Do a randomized hill-climbing search
while (~done && (currentPenalizedScore > maxScore || invalidOperation))
    invalidOperation = false;
    maxScore = currentPenalizedScore;
    % Check all possible remaining candidate edges
    [search, numevals] = search.EvalAllEdges();
    evalcount = evalcount + numevals;
    candidateEdges = find(search.lldiff > bf_thresh);
    % If no more candidate edges, exit
    if (isempty(candidateEdges))
        done = true;
        break;
    end
    
    % Randomly pick a candidate edge
    currentEdge = randsample(candidateEdges, 1);
    [choseni, chosenj] = ind2sub(size(search.llscore), currentEdge);
    newEdgeScore = search.llscore(choseni, chosenj);
    newEdgeDiffScore = search.lldiff(choseni, chosenj);
    oldEdgeScore = newEdgeScore - newEdgeDiffScore;

    % Compute number of free parameters for candidate network structure
    [search, newNumParameters] = search.GetNumParameters(choseni, chosenj, true);
    % Compute BIC penalty score
    newPenaltyScore = (0.5 * newNumParameters) * log(numSample);
    if (~searchParameter.bic)
        % Compute AIC penalty score
        newPenaltyScore = newNumParameters;
    end
    newScore = currentScore - oldEdgeScore + newEdgeScore;
    newPenalizedScore = newScore - newPenaltyScore;
    % Decide if a candidate edge should be added. If so, add edge.
    if (newPenalizedScore > currentPenalizedScore)
        [search, success] = search.DoEdgeAdd(choseni, chosenj, newEdgeDiffScore);
        if (success)
            % Debugging
            if (verbose)
                outline = ['Edge op succesful;', num2str(choseni), ';', num2str(chosenj), ';', num2str(newEdgeDiffScore)];
                disp(outline);
            end
            currentScore = newScore;
            currentPenalizedScore = newPenalizedScore;
            llscores(end+1) = llscores(end) + currentEdgeScore;
            addedeges{end+1} = [choseni, chosenj];
            numedges(end+1) = search.edges;
            evalnums(end+1) = evalcount;
        else
            % This edge operation is invalid: creates a cycle or exceedes number of parents
            invalidOperation = true;
        end
    else
        % Stop search. This edge operation (add, delete, or reverse) did not improve BIC score
        newPenalizedScore - currentPenalizedScore;
    end
end

%% Update edge weight for continuous to continuous nodes to reflect regression coefficient sign 
search = search.UpdateCoefficients();

%% Set up output to be similar to learnHybridBayesNet()
BayesNet.adjMatrix = search.bn;
BayesNet.weightMatrix = search.weightMatrix;

outstats.lldiffs = llscores;
outstats.numedges = numedges;
outstats.numevals = evalnums;
outstats.addededge = addedeges;
outstats.search = search;
