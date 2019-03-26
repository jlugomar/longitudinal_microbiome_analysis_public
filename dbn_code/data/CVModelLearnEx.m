function [auc, numnodes, classacc, cvAdjMat] = CVModelLearnEx(data, cols, pheno, priorPrecision, ...
    folds, experimentname, verbose, randfolds, ALG, bffilter)
%[auc, numnodes, classacc] = CVModelLearn(data, cols, pheno, priorPrecision, ...
%       folds, experimentname, verbose, ALG)
% Does crossvalidation of model building and testing on data.  Will build a
% different network for each fold of cross validation.  Provides a
% good test of various parameter settings if called with differen values of
% PRIORPRECISION.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%       priorPrecision.BFTHRESH: minimum increase in Log Likelihood for 
%           inclusion of an edge in the network.  Deafult = 0;
% FOLDS: Number of folds in the cross-validation to perform.  Default = 5.
% EXPERIMENTNAME: string that will be used in fileoutput names.  Should
%   represent a valid filename
% VERBOSE: boolean.  If true, increases output.
% RANDFOLDS: boolean, if true, will randomly select each cross validation
%   fold.  not gauranteed to be the same size.  (default = true)
% ALG: Network Learning algorithm indicator (optional)
%   ALG == 1 : use K2 (defualt) search for building networks
%   ALG == 2 : use PhenoCentric search for building networks
%   ALG == 3 : use Exhaustive search for building networks
% BFFILTER : if true, will first filter out any SNPs that don't have
%   required BayesFactor Threshold with the phenotype; but does this within
%   the cross-validation loop.  Dafualt = false.
%
% OUTPUT: 
% AUC: the final AUC of the exercise; aggregated over each fold and
%   combined for the testing set of each fold.
% NUMNODES: size of each fold, in number of variables. 
% CLASSAC: accuracy per class of the phenotype.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    folds = 5;
end
if (nargin < 6)
    experimentname = 'bayesnet-CV';
end
if (nargin < 7)
    verbose = true;
end
if (nargin < 8)
    randfolds = true;
end
if (nargin < 9)
    ALG = 1;
end
if (nargin < 10)
    bffilter = false;
end

if (~isfield(priorPrecision,'BFTHRESH'))
    BFTHRESH = 0;
else
    BFTHRESH = priorPrecision.BFTHRESH;
end

% find pheno col; count data
phencol = strmatch(pheno, cols, 'exact');
[ncases,ncols] = size(data);

if (randfolds)
    % randomly split data into N-fold CV sets:
    r = ceil(rand(1,ncases) * folds);
else
    r = zeros(1,ncases);
    for i = 1:folds
        r(ceil(ncases/folds)*(i-1)+1:ceil(ncases/folds)*i) = i;
    end
    r = r(1:ncases);
end
auc = 0;
numnodes = zeros(1,folds);

basecols = cols;
cvPClass = [];
cvPredZs = [];
cvTrueClass = [];
cvAdjMat = zeros(length(cols),length(cols));
% for each fold in teh cross-validation
for k = 1:folds
    cvdata = data(r ~= k,:);
    cvtest = data(r == k,:);
    if (verbose)
        fprintf(1,'Starting Fold %d!\n',k);
    end
    if (isempty(cvtest))
        if (verbose)
            fprintf(1,'Skipping Fold %d because there is no test data!\n',k);
        end
        continue;
    end

    % check values of discrete vars::
    d = IsDiscrete(data, 5);
    discvals = cell(size(d));
    for i = 1:length(cols)
        if (d(i))
            discvals{i} = num2cell(unique(data(:,i),'legacy'));
        else 
            discvals{i} = {};
        end
    end

    % learn a BN for this fold of data:
    if (bffilter)
        [~,MBNet] = BFFilterBNLearn(cvdata, cols, pheno, ALG, BFTHRESH, verbose, priorPrecision);
        % fix the MBNet.mb field to point at indices into actual columns
        % here, which is good for MBNet.AddDiscVals()
        for i = 1:length(MBNet.cols)
            match = strcmp(MBNet.cols{i},cols);
            MBNet.mb(i) = find(match);
        end
    else
        if (ALG == 3)
            cols = basecols;
            BN = FullBNLearn(cvdata, cols, pheno, BFTHRESH, ...
                [experimentname,'-fold',num2str(k),'_exhaustive_1'], priorPrecision, d);
            MBNet = BN.MakeIntoMB();
        elseif (ALG == 2)
            MBNet = LearnPhenoCentric(cvdata, cols, pheno, priorPrecision, BFTHRESH, verbose, d);
        else
            MBNet = LearnStructure(cvdata, cols, pheno, priorPrecision, [experimentname,'-fold',num2str(k)], verbose);
        end
    end
    % sometimes hit a problem here that the FullMB can have no nodes left
    MBNet = MBNet.AddDiscVals(discvals);
    
    numnodes(k) = length(MBNet.mb)-1;
    if (numnodes(k) == 0)
        % no network learned; no associations worth making
        acc = .5;
        p = ones(size(cvtest,1),1);
        z = p * .5;
        if (verbose)
            fprintf(1,'No Network worth making\n');
        end
    else
        % there was a network
        if (verbose)
            fprintf(1,'Learning Network Parameters\n');
        end
        if (~bffilter)
            % Parameters already learned with the BFFilterBNLearn function
            MBNet = LearnParamsBN(MBNet);
        end
        if (verbose)
            fprintf(1,'Predicting on Testing Data\n');
        end
        MBTestNet = MBNet.ReplaceData(cvtest,cols);
%        MBTestNet = MBNet.ReplaceData(cvdata,cols);
        [acc, p, z] = PredictPheno(MBTestNet, verbose);
    end
    FullNet = MBNet.InflateFromMBtoFull(cols);
    cvAdjMat = cvAdjMat + FullNet.adjmat;
    cvPClass = [cvPClass; p];
    cvPredZs = [cvPredZs; z];
    cvTrueClass = [cvTrueClass; cvtest(:,phencol)];
end
cvAdjMat = cvAdjMat ./ folds;
[auc, classacc] = AUCWorker(acc,cvPClass,cvPredZs,cvTrueClass, true, false, verbose);
    
