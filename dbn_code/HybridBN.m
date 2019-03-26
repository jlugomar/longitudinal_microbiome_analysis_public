function convexhullAUC = HybridBN(fname, pheno, analysis_title, delim)
%convexhullAUC = HybridBN(fname, pheno, analysis_title)
% This function is a good starting point for CG Bayesian Network analysis.
% 
% Will build the best bayes net and measure its performance for the file
% given in FNAME.
% This should be a file with columns representing each variable and rows
% representing each data record.
%
% INPUT:
%   FNAME : the file name, a text file tab-delimited that has column
%       names for each variable.
%   PHENO : a string representing the variable to predict.  Must match one
%       of the columns of the file, above
%   ANALYSIS_TITLE : a descriptive title to use for this analysis.  Will
%       output network files with this as the filename root.
%   DELIM : file delimiter string. Common values : '\t' (tab, for text 
%       files) and ',' (comma, for csv files).  Optional. Defaults to any 
%       whitespace character if omitted.
%
% OUTPUT:
%   CONVEXHULLAUC : the AUC of the prediction of the learned network.  The
%       convex hull of the ROC is given as the AUC.  This is a suitable
%       measure of the predictive accuracy of the network.
% 
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


fprintf(1,'\n\n======  start new analysis: %s  =====\n\n',analysis_title);
bnpathscript;
tic;
    

% common parameter values:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 3;

verbose = true;
doAUC = true;

fprintf(1,'Reading data from %s\n', fname);
if (nargin < 4)
    [data, cols] = RCSVLoad(fname, true);
else
    [data, cols] = RCSVLoad(fname, false, delim);
end    
fprintf(1,'Learning Most Predictive Network Structure for %s\n', analysis_title);
MBNet = LearnStructure(data, cols, pheno, priorPrecision, [analysis_title,'-net']);
fprintf(1,'Learning Network Parameters\n');
[tree, nodes] = LearnParams(MBNet, '', data, cols, priorPrecision);
fprintf(1,'Predicting on Training Data\n');
[acc, p, z] = PredictPheno(tree, nodes, '', pheno, data, cols, verbose);

% output AUC or accuracy measure:
phncol = strmatch(pheno, cols, 'exact');
convexhullAUC = AUCWorker(acc,p,z,data(:,phncol),doAUC);


