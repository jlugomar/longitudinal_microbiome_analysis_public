function [auc, good, numvars, testauc, model, classacc, testclassacc] = ...
    LogFitAndTest(data, cols, testdata, testcols, pheno, verbose)
%[auc, good, numvars, testauc, model, classacc, testclassacc] = LogFitAndTest(data,cols, testdata,testcols, pheno,verbose)
%
% Builds a logistic regression model, building on data and testing on testdata.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% TESTDATA: testing data array
% TESTCOLS: testing column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% VERBOSE: boolean.  If true, increases output.
%
% OUTPUT: 
% AUC: the training AUC of the model on the data.
% GOOD: boolean array, parallel to COLS, indicating the variables that are 
%   used in the model.
% NUMVARS: size of model, in number of variables. 
% TESTAUC: AUC of model on test dataset.
% MODEL: Column names of variables used in the model.
% CLASSACC: accuracy of model on training data, per class of the phenotype.
% TESTCLASSACC: accuracy of model on test data, per class of the phenotype.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


TRI_TEST_STRAT = 0;

if (nargin < 6)
    verbose = true;
end

% find pheno col
phencol = strmatch(pheno, cols, 'exact');
xind = true(size(cols));
xind(phencol) = false;
y = data(:,phencol);
x = data(:,xind);

% find pheno col in test data
testphencol = strmatch(pheno, testcols, 'exact');
testxind = true(size(testcols));
testxind(testphencol) = false;
testy = testdata(:,testphencol);
testx = testdata(:,testxind);

% trinary test strategy 
if (TRI_TEST_STRAT == 1)
    % in training set kill cols with class == 0;
    % in testing set re-assign cols with (class == 0) to (class = 1).
    x = x(y ~= 0,:);
    y = y(y ~= 0);
    y(y == -1) = 0;
    testy(testy == 0) = 1;
    testy(testy == -1) = 0;
end


% fit logistic regression w/no 2nd-order interaction terms:
if (verbose)
    fprintf(1,'Learning Logistic Model:\n');
end
[b,dev,stats] = glmfit(x,y,'binomial','link','logit');

% use the values that are statistically significant:
good = stats.p(2:end) <= 0.05;
bgood = stats.p <= 0.05;
numvars = sum(good);
model = cols(good);
if (stats.p(1) <= 0.05)
    const = 'on';
else
    const = 'off';
end

% predict on the data
if (verbose)
    fprintf(1,'Computing Logistic Predictions on Training Data:\n');
end
yzs = glmval(b(bgood),x(:,good),'logit','constant',const);
testyzs = glmval(b(bgood),testx(:,good),'logit','constant',const);
if (size(yzs,2) > 1)
    yzs = yzs(:,2);
end
if (size(testyzs,2) > 1)
    testyzs = testyzs(:,2);
end
ypred = yzs > .5;
acc = ypred == y;
testypred = testyzs > .5;
testacc = testypred == testy;
% compute AUC of these predictions:
[auc,classacc] = AUCWorker(acc,ypred,yzs,y,true,true,verbose);
[testauc,testclassacc] = AUCWorker(testacc,testypred,testyzs,testy,true,true,verbose);
