function DBN = learnDynamicBayesNetwork(filename, pheno, cols2learn, maxParents, localNetStruc, intraEdges, mle, nu, alpha, sigma2, BFTHRESH, MAXDISCVALS)
%DBN = learnDynamicBayesNetwork(filename, pheno, maxParents, localNetStruc, intraEdges, mle, nu, alpha, sigma2, BFTHRESH, MAXDISCVALS)
%
% Constructs a Dynamic Bayesian Network (DBN). 
% Outputs to both trivial graph format and SIF format for Cytoscape.
%
% INPUT:
%     filename      % filename with data for learning a DBN
%     pheno         % string matching one of the colunms. Defaults to none. (optional)
%     cols2learn    % variables to be learned. Input as indices. Defaults to all. (optional)
%     maxParents;   % maximum number of parents for any node. Defaults to 3 max parents per node. (optional)
%     localNetStruc % if true, compare only local network structures using log Bayes factor instead of BIC/AIC score. Defaults to false. (optional)
%     intraEdges;   % if true, intra edges will be learned in the network structure. Defaults to false. (optional)
%     mle;          % if true, maximum likelihood estimation will be used. 
%                   % Otherwise a Bayesian estimation (MAP) will be employed. Defaults to true. (optional)
%     nu            % prior sample size for prior variance estimate. Defaults to 1. (optional)
%     alpha         % prior sample size for discrete nodes. Defaults to 1. (optional)
%     sigma2        % prior variance estimate. Defaults to 1. (optional)
%     BFTHRESH      % Bayes factor threshold for local network structure learning algorithm. Dafaults to 0. (optional)
%     MAXDISCVALS   % maximum number of unique values a column can contain before it can no longer be considered discrete. Appropriate values are
%                   % dependent on the sample size, but probably between 3 and 8. Defaults to 4; (optional)
%
% OUTPUT:
%     DBN           % BayesNet class object representing a 2-stage dynamic bayes net (TSBN)
%
% Sample calls:
% G = learnDynamicBayesNetwork('example', '', [], 3, false, false, true); %MLE w/o intra edges
% G = learnDynamicBayesNetwork('example', '', [], 3, false, true, true); %MLE w/ intra edges
% G = learnDynamicBayesNetwork('example', '', [], 3, true, false, false, 1, 1, 1); %MAP w/o intra edges
% G = learnDynamicBayesNetwork('example', '', [], 3, true, true, false, 1, 1, 1); %MAP w/ intra edges
%
    if (nargin < 12)
        MAXDISCVALS = 4;
    end
    if (nargin < 11)
        BFTHRESH = 0;
    end
    if (nargin < 10)
        sigma2 = 1;
    end
    if (nargin < 9)
        alpha = 1;
    end
    if (nargin < 8)
        nu = 1;
    end
    if (nargin < 7)
        mle = true;
    end
    if (nargin < 6)
        intraEdges = false;
    end    
    if (nargin < 5)
        localNetStruc = false;
    end
    if (nargin < 4)
        maxParents = 3;
    end
    if (nargin < 3)
        cols2learn = [];
    end
    if (nargin < 2)
        pheno = '';
    end
    if (nargin < 1)
        filename = 'DBN_test';
    end
    
    %% Load data
    fname = [filename, '.tsv'];
    [numdata, cols] = RCSVLoad(fname, false, '\t');
    data = numdata(:, 2:end);
    cols = cols(2:end);
    subjectIDs = numdata(:, 1);
    numVariables = size(data, 2);
    
    %% Take data and split into two (overlapping) time points
    [d0, dn, nsids] = MakeTSBNData(data, subjectIDs);

    %% Determine which variables are to be treated as discrete
    disc = IsDiscrete(data, MAXDISCVALS);

    %% Determine the set variables to be learned
    if isempty(cols2learn)
        cols2learn = 1:numVariables;
    end
    targetCols = false(1, numVariables);
    targetCols(cols2learn) = 1;

    %% Make main call to learn structure of the TSBN dataset
    priorPrecision.nu = nu;
    priorPrecision.alpha = alpha;
    priorPrecision.sigma2 = sigma2;
    priorPrecision.maxParents = maxParents;
    priorPrecision.mle = mle;
    searchParameter.d0 = d0;
    searchParameter.dn = dn;
    searchParameter.DBN = true;
    searchParameter.nophenotype = true;
    searchParameter.annealing = false;
    searchParameter.local = localNetStruc;
    searchParameter.bic = mle;
    searchParameter.unwrapped = false;

    if (~intraEdges)
        % FullBNLearn() will write TSBN to GML and TGF if verbose is TRUE
        [BNet, outstats] = FullBNLearn(data, cols, pheno, BFTHRESH, filename, ...
                           priorPrecision, disc, false, searchParameter, targetCols);
        
        % Convert TSBN into DBN:
        newdata = [d0, dn];
        DBN = BNet.ConvertToDBN(newdata);

        % Adding suffix _dbn to properly label DBN
        dbn_filename = [filename, '_dbn'];
        
        % Write DBN to GML and TGF
        DBN.WriteToGML(dbn_filename);
%        DBN.WriteToTGF(dbn_filename);
    else
        % Adding suffix _dbn to properly label DBN
        dbn_filename = [filename, '_dbnIntra']; 

        searchParameter.unwrapped = true;
        % FullBNLearn() will write DBN directly to GML and TGF if verbose is TRUE
        [DBN, outstats] = FullBNLearn(data, cols, pheno, BFTHRESH, dbn_filename, ...
                          priorPrecision, disc, true, searchParameter, targetCols);
    end
