function [convexhullAUC, classAcc, varargout] = AUCWorker(acc,p,z,ptrue,doAUC,dontInvert, verbose, doplot)
%[convexhullAUC, classAcc, varargout] = AUCWorker(acc,p,z,ptrue,doAUC,dontInvert, verbose, doplot)
%
% Main function for computing AUCs from predictions on data.  Works for
% Bayesian Network classifiers and others provided in the CGBayesNets
% pacakge (Niave Bayes, Logistic Regression).  Will provide confidence
% intervals and draw AUC plots.
%
% INPUT:
% ACC : accuracy in raw % correct
% P : class predictions of the classifier.  Parallel to Z and PTRUE.
% Z : confidence of class prediction, in [0,1]. Parallel to P and PTRUE.
% PTRUE : true class. Parallel to Z and P.
% DOAUC : boolean, if true, use full AUC.
% DONTINVERT: true if used with some other prediction mechanism than the
%   Chang-McGeachie BayesNet prediction.  The BN needs to have its
%   predictions inverted on negative examples.  Other stuff (logistic
%   regression models) does not.
% VERBOSE: if true, writes out more info to the screen.
% DOPLOT: if true, will pop up a figure of the AUC and ConvexHull AUC.
%
% OUTPUT:
% CONVEXHULLAUC: The convex-hull AUC of the classifier
% CLASSACC: The accuracy (in % correct) broken up per class of the phenotype.
%   Useful for determining if the accuracy comes from predicting on class
%   well at the expense of the other(s).
% VARARGOUT{1}: rawAUCConf: the confidence interval of the raw AUC value.
% VARARGOUT{2}: convexhullAUCConf: the confidence interval of the convex
%   hull AUC.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 8)
    doplot = false;
end
if (nargin < 7)
    verbose = true;
end
if (nargin < 6)
    dontInvert = false;
end
if (nargin < 5)
    doAUC = true;
end

if ((length(unique(ptrue,'legacy')) > 2) || (length(unique(ptrue,'legacy')) == 1))
    doAUC = false;
end
pclasses = unique(ptrue, 'legacy');
pclasses = sort(pclasses);
if (length(pclasses) == 2 && pclasses(1) ~= 0 && pclasses(2) ~= 1)
    % switch p and ptrue to be in terms of 0 and 1:
    ptemp = p;
    ptemp(p == pclasses(1)) = 0;
    ptemp(p == pclasses(2)) = 1;
    p = ptemp;
    ptemp = ptrue;
    ptemp(ptrue == pclasses(1)) = 0;
    ptemp(ptrue == pclasses(2)) = 1;
    ptrue = ptemp;
end

classAcc.classes = unique(ptrue,'legacy');
classAcc.acc = zeros(size(classAcc.classes));
for i = 1:length(classAcc.classes)
    c = classAcc.classes(i);
    inds = ptrue == c;
    if (~isempty(c))
        classAcc.acc(i) = sum(p(inds) == ptrue(inds))/length(ptrue(inds));
    end
end

if (doAUC)
    if (~dontInvert)
        % first invert the predictions on the first class:
        % (only makes sense for binary phenotype)
        z(p == 0) = 1 - z(p == 0);
    end
    % then call the ROC routines to compute true positives, false positive rates:
    [tp,fp] = roc(ptrue, z);
    if (doplot)
        rocplot(tp,fp);
    end
    rawAUC = auroc(tp,fp);
    % and convex-hull of ROC:
    [tp,fp] = rocch(ptrue, z);
    convexhullAUC = auroc(tp,fp);
    % variance of AUC
    var = aucvar(ptrue,z);
    std = sqrt(var);
    % get conf intervals
    conf = 0.95;
    rawAUCConf = norminv([(1-conf)/2, conf + (1-conf)/2], rawAUC, std);
    convexhullAUCConf = norminv([(1-conf)/2, conf + (1-conf)/2], convexhullAUC, std);

    if (verbose)
        fprintf(1,['\tPrediction AUC: ',num2str(rawAUC),' [', num2str(rawAUCConf(1)),',',...
            num2str(rawAUCConf(2)),'], and convex-hull AUC: ', num2str(convexhullAUC),...
            ' [', num2str(convexhullAUCConf(1)),',', num2str(convexhullAUCConf(2)),'].\n']);
    end
    
    if (nargout > 2)
        varargout{1} = rawAUCConf;
        varargout{2} = convexhullAUCConf;
    end
else
    acc = sum(acc)/length(acc);
    if (verbose)
        fprintf(1,['\tPrediction Accuracy: ', num2str(acc), '\n']);
    end
    convexhullAUC = acc;
end

