function [BayesNet, outstats] = SimAnnealLearn(contData, discData, priorPrecision, phencol, searchParameter, search)
%[BayesNet, outstats] = SimAnnealLearn(contData, discData, priorPrecision, phencol, searchParameter, search)
%
% Performs a simulated annealing search over the entire network, considering
% one move at a time, which is either adding or removing an randomly chosen
% edge.  Move is accepted with probability 1 if it improves the loglikelihood
% and with probability ~ exp(Delta-SEARCHPARAMETER.BF_THRESH)/T); where
% Delta is the change in loglikelihood of the move; where T is = 
% (N^3 - count) / N^3; where N is the number of nodes in the data; and count 
% is the number of moves accepted so far.
% Roughly O(n^3).  
%
% Input:
% CONTDATA: numeric data matrix of continous nodes (indexed by node,sample)
% DISCDATA: numeric data matrix of discrete nodes (indexed by node,sample)
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% PHENCOL: column index for the phenotype column.
% SEARCHPARAMETER: structure with the following fields:
%   searchParameter.BF_THRESH: a Bayes Factor threshold for including an 
%       edge.  Edges with likelihood below BF_THRESH will not be included.  
%       This is the main stopping criterion.
%   searchParameter.nophenotype: a boolean, if true, means that the
%       network should be built without any particular phenotype
%   searchParameter.SA_Temp_Mult: if present, will use this as a 
%       multiplier of the usual simulated annealing temperature schedule, 
%       making the search take (SA_Temp_Mult) times as long as normal.
%   searchParameter.Back_Step_Mult: if present, will use this as a 
%       multiplier of the probability of making a backwards move. 
%   searchParameter.DBN: a boolean, if true, indicates that we're learing a
%       two-stage dynamic bayes net, and will use the TSDBNStateTrackerSearch
%       class. Default = false.
%       A dynamic bayes net must further define:
%   searchParameter.t0cont: continuous vars of the first of two-stages
%   searchParameter.tncont: continuous vars of the second of two-stages
%   searchParameter.t0disc: discrete vars of the first of two-stages
%   searchParameter.tndisc: discrete vars of the second of two-stages        
% 
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
%   outstats.temp: simulated annealing temperature at each step in the
%       network search.
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

% consider using a multiplier of acceptance probability:
if (~isfield(searchParameter,'Back_Step_Mult'))
    BACK_STEP_MULT = 10;
else
    BACK_STEP_MULT = searchParameter.Back_Step_Mult;
end
% multiples of n^3 for max number of steps to take: 
if (~isfield(searchParameter,'SA_Temp_Mult'))
    SA_TEMP_MULT = 0.2;
else
    SA_TEMP_MULT = searchParameter.SA_Temp_Mult;
end

% phenotype node can't have any edges going in:
if (~isfield(searchParameter,'nophenotype'))
    nophenotype = false;
else
    nophenotype = searchParameter.nophenotype;
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

bf_thresh = searchParameter.BF_THRESH;
backtracking = true;
checkRepeats = false;

% set up the search:
if (nargin < 6)
    if (DBN)
        if (unwrapped)
            search = UnwrappedTSDBNSTS(contData, discData, priorPrecision, phencol, ...
                backtracking, bf_thresh, nophenotype, checkRepeats);
            search = search.Init();
            search = search.SetUnwrappedTimeSeries();
        else
            % a Dynamic BN must define:
            % searchParameter.t0cont
            % searchParameter.tncont
            % searchParameter.t0disc
            % searchParameter.tndisc
            % 
            % set up the search to do a dynamic bayesian network:
            search = TSDBNStateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats);
            search = search.Init();
            search = search.SetTimeSeries(searchParameter.t0cont, searchParameter.t0disc, ...
                searchParameter.tncont, searchParameter.tndisc);
        end
    else
        search = StateTrackerSearch(contData, discData, priorPrecision, phencol, ...
            backtracking, bf_thresh, nophenotype, checkRepeats);
        search = search.Init();
    end
else
    search.contData = contData;
    search.discData = discData;
    search.priorPrecision = priorPrecision;
    search.backtracking = backtracking;
    search.bf_thresh = bf_thresh;
    search.nophenotype = nophenotype;
    search.checkRepeats = checkRepeats;
    search = search.ReOpenSearch();
    if (DBN && ~unwrapped)
        % this should work with TSDBNStateTrackerSearch class DBNs
        search.t0cont = searchParameter.t0cont;
        search.tncont = searchParameter.tncont;
        search.t0disc = searchParameter.t0disc;
        search.tndisc = searchParameter.tndisc;
%    elseif (DBN && unwrapped)
        % don't set up those parameters for an unwrapped DBN
        
    end
end

lldiffs = 0;
numedges = 0;
addedeges = {};
temp = 0;
done = false;
count = 0;
evalcount = 0;
evalnums = 0;
%% do selection loop:
while (~done && ~search.done)
    count = count + 1;
    
    % get list of possible moves (edges to add/remove):
    accept = false;
    SA_Temp = BACK_STEP_MULT * max(0, SA_TEMP_MULT * length(search.bn)^3 - count) ...
        / (SA_TEMP_MULT * length(search.bn)^3);
    if (SA_Temp == 0)
        % once the temp hits zero, we're done searching around and will
        % stop at any previously-visited state.
        search.checkRepeats = true;
        % check for final doneness: if all edges are negative and SA_Temp == 0.
        [search, numevals] = search.EvalAllEdges();
        evalcount = evalcount + numevals;
        possible_edges = find(search.lldiff > bf_thresh);
    else
        % possible edges are to add a new edge that isn't closed off;
        % or to remove an existing edge from the bn:
        possible_edges = find(~search.closed | search.bn);
    end
    % if nothing to do, quit:
    if (isempty(possible_edges))
        done = true;
        break;
    end
    
    % randomly check each of these moves until we find one that we accept
    porder = randperm(length(possible_edges));
    i = 0;
    while (~accept && i < length(porder))
        i = i + 1;
        pe = possible_edges(porder(i));
        [choseni,chosenj] = ind2sub(size(search.lldiff), pe);
        [search, evalq] = search.EvalOneEdge(choseni, chosenj);
        if (evalq)
            evalcount = evalcount + 1;
        end
        lldiffchosen = search.lldiff(choseni, chosenj);
        accept = rand(1) < exp((lldiffchosen-bf_thresh)/SA_Temp);
    end
    
    if (~accept && (SA_Temp == 0))
        done = true;
    end
    
    if (~done && accept)
        [search, success] = search.DoEdgeAdd(choseni, chosenj, lldiffchosen);
        if (success)
            lldiffs(end+1) = lldiffs(end) + lldiffchosen;
            addedeges{end+1} = [choseni,chosenj];
            numedges(end+1) = search.edges;
            temp(end+1) = SA_Temp;
            evalnums(end+1) = evalcount;
        end
    end
end


%% set up output to be similar to learnHybridBayesNet()
BayesNet.adjMatrix = search.bn;
BayesNet.weightMatrix = search.weightMatrix;

outstats.lldiffs = lldiffs;
outstats.numedges = numedges;
outstats.numevals = evalnums;
outstats.temp = temp;
outstats.addededge = addedeges;


