function convexhullAUC = DoTestNetworkEval(nodesketch, datafile, pheno, ...
    priorPrecision, networkname, doAUC, verbose)
% convexhullAUC = DoTestNetworkEval(nodesketch, datafile, pheno, priorPrecision, networkname, doAUC, verbose)
% 
% DEPRICATED
%
% Takes a nodesketch representing a BN and makes calls to the
% functions to learn parameters and perform testing.  Doesn't learn the
% network structure.
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 7)
    verbose = true;
end
if (nargin < 6)
    doAUC = true;
end

fprintf(1,'Reading data from %s\n', networkname);
[data,cols] = ReadHeaderDataFile(datafile);
fprintf(1,'Learning Network Parameters\n');
[tree, nodes] = LearnParams(nodesketch, '', data, cols, priorPrecision);
fprintf(1,'Predicting on Training Data\n');
[acc, p, z] = PredictPheno(tree, nodes, '', pheno, data, cols, verbose);

% output AUC or accuracy measure:
phncol = strmatch(pheno, cols, 'exact');
convexhullAUC = AUCWorker(acc,p,z,data(:,phncol),doAUC);

