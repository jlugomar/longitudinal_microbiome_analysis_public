function [MBnet, FullBN, outstats] = LearnStructure(data, cols, pheno, priorPrecision, outfilename, ...
    verbose, disc)
%[MBNet, FullBN, outstats] = LearnStructure(data, cols, pheno, priorPrecision, outfilename, verbose, disc)
%
% Main function for learning the structure of a Bayes Network; following
% the K2 strategy of only making parents between nodes with lower Bayes 
% Factor than the current node.  This can be called directly or invoked with
% BFFilterBNLearn().
%
% INPUT:
%   DATA: data array (data points in rows by variables in columns)
%   COLS: column names, a cell array of strings
%   PHENO: a string representing the phenotype column to predict.  Is matched
%     against the COLS array.
%   PRIORPRECISION : a structure indicating the several parameter settings
%     for the prior of the bayesian search.  These include:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%       priorPrecision.BFTHRESH: minimum increase in Log Likelihood for 
%           inclusion of an edge in the network.  Deafult = 0;
%   OUTFILENAME : filename for printing out the network structure in a .tgf
%     file.  (Trivial Graph Format).
%   VERBOSE : increases output printed to stdout. (optional, default =
%     true)
%   DISC : (optional) logical array parallel to cols specifying which are
%       discrete.
%
% OUTPUT:
% MBnet: Class BAYESNET object representing the Markov Blanket of the
%   Bayes Net.  Reorders the COLS for this MB.
% FullBN: Class BAYESNET object representing the full network determined 
%   by the K2 algorithm.
% OUTSTATS: structure with fields describing the characteristics of the
%       search procedure, in arrays per "step;" a step is either an edge
%       added or removed.
%   outstats.lldiffs: difference in loglikelihood at each step of the algorithm.
%   outstats.numedges: number of edges in the network at each step of the algorithm.
%   outstats.numevals: number of potential network states evaluated at each
%       step.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    verbose = true;
end
if (~isfield(priorPrecision, 'maxParents'))
    priorPrecision.maxParents = 3;
end
rootNodeName = pheno;


%% figure out which columns (variables) are discrete
if (nargin < 7)
    disc = IsDiscrete(data);
end

origcols = cols;
%% find phenotype column:
inds = 1:length(cols);
phninds = strcmp(pheno, cols);
phncol = inds(phninds);
% move phenotype coln to the last position:
moverindex = true(size(cols));
moverindex(phncol) = false;
data = [data(:,moverindex), data(:,phncol)];
cols = cols(moverindex);
cols{end+1} = pheno;
phendisc = disc(~moverindex);
disc = disc(moverindex);
disc(end+1) = phendisc;


%% reshape the dataformat for network learning
% according to conventions of learning; we should transpose the data here
% discData: rows are variables, columns are sample observations 
discData = single(data(:,disc)'); 
discNodeNames = cols(disc);
% flip continuous cols, too
contData = data(:,~disc)';
contNodeNames = cols(~disc);
numContNodes = size(contData,1);

% other parameters for learning bayesnets:
initialBayesNet = []; % initial Bayes network     
searchParameter.STEPWISE=true;
searchParameter.MAX_NUM_PARENTS=priorPrecision.maxParents;
if (isfield(priorPrecision,'BFTHRESH'))
    searchParameter.BFTHRESH=priorPrecision.BFTHRESH;
end



%% filter variables 
% sort the nodes according to the likelihood of being independent of the rest nodes
if (verbose)
    fprintf(1,'Start sorting the nodes likelihood of independence from other nodes\n');
end
tempLogLLH_Ind = zeros(1,numContNodes);
for v=1:numContNodes
    [tempModel]=learnLocalBN_DiscToCont([], contData(v,:), priorPrecision);    
    tempLogLLH_Ind(v) = tempModel.logLLH;
end
numevals = numContNodes;
[tempLogLLH_Ind,tempIdx]=sort(tempLogLLH_Ind,'ascend');
contData = contData(tempIdx,:);
contNodeNames = contNodeNames(tempIdx);
if (verbose)
    fprintf(1,'!Done filtering and reordering nodes by independence from other nodes!\n');
end

%% learn Bayes network
if (verbose)
    fprintf(1,'Start learning Bayesian network structure!\n');
end
% if we have a continuous phenotype, include its index relative to the
% contData array
if (~phendisc)
    % find phencol in the contNodeNames list
    contphenind = find(strcmp(contNodeNames, pheno));
else
    contphenind = [];
end

[BN,outstats]=learnHybridBayesNet(contData,discData,priorPrecision,contphenind, ...
    searchParameter,initialBayesNet);
outstats.numevals = outstats.numevals + numevals;
%note: BN.adjMatrix(m,n)=1 if m is a parent of n
if (size(contNodeNames,2) > size(contNodeNames,1))
    if (size(discNodeNames,2) > size(discNodeNames,1))
        nodeNames = [contNodeNames';discNodeNames'];
    else
        nodeNames = [contNodeNames';discNodeNames];
    end
else
    if (size(discNodeNames,2) > size(discNodeNames,1))
        nodeNames = [contNodeNames;discNodeNames'];
    else
        nodeNames = [contNodeNames;discNodeNames];
    end
end
%note: after Bayes net learning, the whole node names start with contNodes followed by discNodes
if (verbose)
    fprintf(1,'!Done learning Bayesian network!\n');
end

%% now make Actual BayesNet structure:
% full bayes net isn't actually worth much in this context, since we only
% computed the markov blanket bayes net.
discmap = [false(1,numContNodes),true(1,length(discNodeNames))];
FullBN = BayesNet([],[outfilename, num2str(priorPrecision.nu)],BN.adjMatrix,...
    BN.weightMatrix,[],false,discmap,[contData;discData]',nodeNames',rootNodeName,priorPrecision,{});
% reconstruct the data columns so they match teh nodeNames:
FullBN = FullBN.ReorderByNames(origcols);
MBnet = FullBN.MakeIntoMB();
if (verbose)
    MBnet.WriteToTGF();
end


