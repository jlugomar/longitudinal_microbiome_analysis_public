function [auc,MBNet] = BFFilterBNLearn(data, cols, pheno, ...
    algorithm, BF_THRESH, verbose, priorPrecision)
%[auc,MBNet] = BFFilterBNLearn(data, cols, pheno, algorithm, BF_THRESH, verbose, priorPrecision)
%
% Main function for building a CG BayesNet.  This first filters the variables
% by using a Bayes Factor threshhold, removing any that do not meet this level
% of dependence upon the phenotype.  This will then learn a
% network using a K2-style algorithm based on the BF of the variables 
% that maximizes Bayesian posterior likelihood of the data given the
% network.  This learns the structure, parameters, and
% then performs inference on the DATA to test the fit of the network;
% producing AUC as output.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array.
% ALGORITHM: parameter that determines the type of search algorithm used to
%       find the best bayesian network
%   ALGORITHM = 1 : markov blanket network search based on K2-style
%       ordering with backtracking.  Default.
%   ALGORITHM = 2 : markov blanket network search based on exhaustive
%       search for combinations of the markov blanket.  Phenocentric.
%   ALGORITHM = 3 : full network search based on greedy search for best
%       edge at each step.
% BFTHRESH: log Bayes Factor (BF) threshold for variable inclusion.  Variables 
%   with BF lower than this will not be part of the network.  Default = 0.
% VERBOSE: if true, increases diagnostic output
% PRIORPRECISION: a structure including the usual CGBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%
% OUTPUT:
% AUC: AUC of the network, trained and predicted on the input DATA.
% MBNET: BayesNet class object representing the CG Bayes Net.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.
    
if (nargin < 4)
    algorithm = 1;
end
if (nargin < 5)
    BF_THRESH = 0;
end
if (nargin < 6)
    verbose = false;
end
if (nargin < 7)
    % parameters for Bayes network learning
    priorPrecision.alpha = 10; %prior frquency for discrete data
    priorPrecision.nu = 10; %prior sample size for continuous data
    priorPrecision.sigma2 = 1; %variance in the prior (hypothetical) samples for continuous data
    priorPrecision.maxParents = 3;
end

phencol = find(strcmp(pheno, cols));

if (verbose)
    fprintf(1,'Computing Bayes Factors\n');
end
BFs = BayesFactorScore(data, cols, pheno, priorPrecision);
keep = (BFs > BF_THRESH);
keep(phencol) = true;
phndata = data(:,phencol);
data = data(:,keep);
cols = cols(keep);
if (verbose)
    fprintf(1,'Found %d variables with BF > %f\n', sum(keep), BF_THRESH);
end

% use one of 3 network learning algorithms:
if (algorithm == 1)
    if (verbose)
        fprintf(1,'Learning Markov Blanket Bayesian Network for Prediction of %s with K2\n', pheno);
    end
    BNet = LearnStructure(data, cols, pheno, priorPrecision, 'K2_mbnet', verbose);
end
if (algorithm == 2)
    if (verbose)
        fprintf(1,'Learning Pheno-centric Bayesian Network\n');
    end
    BNet = LearnPhenoCentric(data, cols, pheno, priorPrecision, BF_THRESH, verbose);
end
if (algorithm == 3)
    if (verbose)
        fprintf(1,'Learning Full Bayesian Network w Exhaustive Greedy Search\n');
    end
    BNet = FullBNLearn(data, cols, pheno, BF_THRESH, 'Full_exhaustive', priorPrecision);
end
if (verbose)
    fprintf(1,'Done learning Bayesian network!\n');
end


%% infer hybrid Bayes network
if (verbose)
    fprintf(1,['Start predicting root node (', pheno,') by the learned Bayesian network!\n']);
end
MBNet = BNet.MakeIntoMB();
MBNet = LearnParamsBN(MBNet);
if (verbose)
    fprintf(1,'Done Learning Conditional Distributions!\n');
end
[predAccuracy, predClass, predZs] = PredictPheno(MBNet);
fprintf(1,['\tPrediction accuracy (raw pct correct): ',num2str(predAccuracy),'\n']);
[auc, ~] = AUCWorker(predAccuracy, predClass, predZs, phndata);

