function DBN = bootstrapLearnDynamicBayesNetwork(filename, pheno, nboots, edgeThreshold, cols2learn, maxParents, localNetStruc, intraEdges, mle, nu, alpha, sigma2, BFTHRESH, MAXDISCVALS, verbose)
%DBN = bootstrapLearnDynamicBayesNetwork(filename, pheno, nboots, edgeThreshold, cols2learn, maxParents, localNetStruc, intraEdges, mle, nu, alpha, sigma2, BFTHRESH, MAXDISCVALS, verbose)
%
% Function to use bootstrapping to generate multiple datasets and learn on
% each a hybrid bayesian network.  
% Outputs an adjacency matrix that contains fractional proportions for each
% edge, where the fraction is the proportion of bootstrap realizations that
% contain the edge.
%
% INPUT:
%     filename      % filename with data for learning a DBN
%     pheno         % string matching one of the colunms. Defaults to none. (optional)
%     nboots        % number of bootstrap realizations to perform
%     edgeThreshold % threshold for edge confidence
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
%     verbose       % if true, will output files of each bootstrap network (optional)
%                   % Default = false.
%
% OUTPUT:
%     DBN           % bootstrapped BayesNet class object representing a 2-stage dynamic
%                   % bayes net (TSBN)
%
% Sample calls:
% G = bootstrapLearnDynamicBayesNetwork('example', '', 10, 0.7, [], 3, false, false, true); %MLE w/o intra edges
% G = bootstrapLearnDynamicBayesNetwork('example', '', 10, 0.7, [], 3, false, true, true); %MLE w/ intra edges
% G = bootstrapLearnDynamicBayesNetwork('example', '', 10, 0.7, [], 3, true, false, false, 1, 1, 1); %MAP w/o intra edges
% G = bootstrapLearnDynamicBayesNetwork('example', '', 10, 0.7, [], 3, true, true, false, 1, 1, 1); %MAP w/intra edges

    if (nargin < 15)
        verbose = false;
    end
    if (nargin < 14)
        MAXDISCVALS = 4;
    end
    if (nargin < 13)
        BFTHRESH = 0;
    end
    if (nargin < 12)
        sigma2 = 1;
    end
    if (nargin < 11)
        alpha = 1;
    end
    if (nargin < 10)
        nu = 1;
    end
    if (nargin < 9)
        mle = true;
    end
    if (nargin < 8)
        intraEdges = false;
    end    
    if (nargin < 7)
        localNetStruc = false;
    end
    if (nargin < 6)
        maxParents = 3;
    end
    if (nargin < 5)
        cols2learn = [];
    end
    if (nargin < 4)
        edgeThreshold = 0.5;
    end
    if (nargin < 3)
        nboots = 1000;
    end
    if (nargin < 2)
        pheno = '';
    end
    
    %% Load data
    fname = [filename, '.tsv'];
    [numdata, cols] = RCSVLoad(fname, false, '\t');
    data = numdata(:, 2:end);
    cols = cols(2:end);
    subjectIDs = numdata(:, 1);
    numVariables = size(data, 2);
    
    adjmat = zeros(2*length(cols));
    weightMatrix = zeros(2*length(cols));

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
    
    %% Generate a bootstrap data realizations
    for b = 1:nboots
        fprintf(1, 'Starting Boot Sample No. %d\n',b);
        
        currFilename = [filename, '_boot', num2str(b)];
        subjects = unique(subjectIDs);
        % Create a sample with replacement from all subjectsIDs 
        % Note that we need to sample all timeppoints for a given subjectID
        rindex = randi(size(subjects,1), size(subjectIDs, 1), 1);
        sids = subjects(rindex, :);
        sampleData = [];
        sampleSubjectIDs = [];
        for i = 1:length(sids)
            si = subjectIDs == sids(i);
            ds = data(si, :);
            sidarray = sids(i) * ones(size(ds(:, 1)));
            sidarray = sidarray + i*10000; %Making sure every iteration of the loop yield different ids
            sampleData = [sampleData; ds];
            sampleSubjectIDs = [sampleSubjectIDs; sidarray];
        end
        
        %% Take sample data and split into two (overlapping) time points
        [d0, dn, nsids] = MakeTSBNData(sampleData, sampleSubjectIDs);
        searchParameter.d0 = d0;
        searchParameter.dn = dn;
        searchParameter.DBN = true;
        searchParameter.nophenotype = true;
        searchParameter.annealing = false;
        searchParameter.local = localNetStruc;
        searchParameter.bic = mle;
        searchParameter.unwrapped = false;
        
        if (~intraEdges)
            % FullBNLearn() will write current TSBN to GML and TGF if verbose is TRUE
            [BNet, outstats] = FullBNLearn(sampleData, cols, pheno, BFTHRESH, filename, ...
                               priorPrecision, disc, verbose, searchParameter, targetCols);

            % Convert TSBN into DBN:
            newdata = [d0, dn];
            DBN = BNet.ConvertToDBN(newdata);
            
            % Adding suffix _dbn to properly label DBN
            dbn_filename = [currFilename, '_dbn'];

            % Write current DBN to GML and TGF
            if (verbose)
                DBN.WriteToGML(dbn_filename);
%                DBN.WriteToTGF(dbn_filename);
            end
        else
            % Adding suffix _dbn to properly label DBN
            dbn_filename = [currFilename, '_dbnIntra'];

            searchParameter.unwrapped = true;
            % FullBNLearn() will write DBN directly to GML and TGF if verbose is TRUE
            [DBN, outstats] = FullBNLearn(sampleData, cols, pheno, BFTHRESH, dbn_filename, ...
                              priorPrecision, disc, verbose, searchParameter, targetCols);  
        end

        % Update global adjacency matrix
        adjmat = adjmat + DBN.adjmat;
        % Update global weight matrix
        weightMatrix = weightMatrix + DBN.weightMatrix;
    end

    %% Update final boostrapped DBN parameters
    adjmat = adjmat ./ nboots;
    bootScore = adjmat;
    weightMatrix = weightMatrix ./ nboots;
    weightMatrix(adjmat < edgeThreshold) = 0;
    adjmat(adjmat < edgeThreshold) = 0;
    adjmat(adjmat >= edgeThreshold) = 1;
    DBN.adjmat = adjmat;
    DBN.weightMatrix = weightMatrix;
    DBN.bootStrapMatrix = bootScore; %To output the bootstrap score instead of the weightmatrix

    %% Adding suffix _dbn to properly label bootstrapped DBN
    dbn_filename = [filename, '_dbnIntraBoot'];
    if (~intraEdges)
        dbn_filename = [filename, '_dbnBoot'];
    end
    
    %% Write final bootstrapped DBN to GML and TGF
    DBN.WriteToGML(dbn_filename);
%    DBN.WriteToTGF(dbn_filename);
