function [predictions, truelikelihood] = inferenceDynamicBayesNetwork(DBN, testdata_filename, cols2predict, outfilename, verbose)
%
% This function takes a 2-stage dynamic BayesNet class object (DBN) with a 
% determined structure. Next, it learns the parameters for DBN using either 
% the data contained in the DBN object or on different data, and then tests 
% the DBN on a subset of variables given the evidence on the other variables.  
%
% INPUT:
%     DBN                   % BayesNet class object representing 2-stage dynamic bayes net
%                             Structure of DBN must have been previously learned
%     testdata_filename     % filename with test data for inference (optional)
%                             Should be in same format as DBN.data and DBN.cols 
%                             If empty, will test on DBN.data
%     cols2predict          % Variables to be predicted. Input as indices into DBN.cols.
%     outfilename;          % filename where predictions and likelihood will be written   
%
% OUTPUT:
%     predictions           % predicted values for each variable in cols2predict
%     truelikelihood        % log probability of observing a less likely outcome for
%                             each variable than the value assigned to that var by the input data.
%
% Sample call: 
%   [p, t] = inferenceDynamicBayesNetwork(G, 'human_vaginal_microbiota_dbn_sample_noalignment_nosplines', [3:25], 'vaginal_testset');

    if (nargin < 5)
        verbose = false;
    end
    if (nargin < 4)
        outfilename = 'DBN_prediction.csv';
    end
    if (nargin < 3)
        cols2predict = [3:length(DBN.cols)/2];
    end
    if (nargin < 2)
        SEPARATE_TEST = false;
        predictions = -1 * ones(size(DBN.data,1), length(DBN.cols)/2);
    else
        SEPARATE_TEST = true;
        %% load test data
        fname = [testdata_filename, '.tsv'];
        [numtestdata, testcols] = RCSVLoad(fname,false,'\t');
        testdata = numtestdata(:,2:end);
        testcols = testcols(2:end);
        subjids = numtestdata(:,1);
        %% take test data and split into two (overlapping) time points
        [td0, tdn, nsids] = MakeTSBNData(testdata, subjids);
        newtestdata = [td0, tdn];
        % this assumes that testcols is in the same order as train data
        newtestcols = DBN.cols;
        predictions = -1 * ones(size(newtestdata,1), length(newtestcols)/2);
    end

%    fsDBNtree = ForwardSimLearnParams(DBN);
%    if (SEPARATE_TEST)
%            DBN = DBN.ReplaceData(newtestdata, newtestcols);
%    end
%    [preds, likelihood, trueppprobs] = ForwardSim(DBN, fsDBNtree);
%    dlmwrite([outfilename, '_predictions.csv'], preds);
%    dlmwrite([outfilename, '_likelihood.csv'], trueppprobs);
    
    z = predictions;
    logevprobs = predictions;
    truelikelihood = predictions;
    verbose = true;
    %% get prediction and true likelihoods using PredictPartial()
    for i = 1:length(cols2predict)
        % identify current target variable (i.e., pheno)
        offset = length(DBN.cols)/2 + cols2predict(i);
        pheno = DBN.cols{offset};
        
        % get Markov blanket for variable to be predicted
        DBN_MB = DBN.MakeIntoMB(pheno);
        
        % learn DBN parameters
        DBN_MB = DBN_MB.LearnParams();

        % get evidence on many of the variable; predict the other variable, 
        % or output their probabilities.
        match = strcmp(pheno, DBN_MB.cols);
        predvars = find(match);
        evvars = 1:length(DBN_MB.cols);
        evvars = evvars(~match);
        if (SEPARATE_TEST)
            DBN_MBTestNet = DBN_MB.ReplaceData(newtestdata, newtestcols);
            if (verbose)
                fprintf(1,'Predicting %s on test data with a Markov Neighborhood of size %d\n', pheno, length(DBN_MBTestNet.cols));
            end
            [predictions(:,cols2predict(i)), z(:,cols2predict(i)),logevprobs(:,cols2predict(i)),truelikelihood(:,cols2predict(i))] = PredictPartial(DBN_MBTestNet, evvars, predvars);
        else
            if (verbose)
                fprintf(1,'Predicting %s on training data with a Markov Neighborhood of size %d\n', pheno, length(DBN_MB.cols));
            end
            [predictions(:,cols2predict(i)), z(:,cols2predict(i)),logevprobs(:,cols2predict(i)),truelikelihood(:,cols2predict(i))] = PredictPartial(DBN_MB, evvars, predvars);
        end
    end

    save([outfilename, '.mat'], 'predictions');
    dlmwrite([outfilename, '_predictions_mb.csv'], predictions);
    dlmwrite([outfilename, '_likelihood_mb.csv'], truelikelihood);
