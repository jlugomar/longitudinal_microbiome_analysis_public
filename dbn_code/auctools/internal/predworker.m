function [tpval] = predworker(tree, nodes, data, cols, pheno, numsims)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

phencol = strcmp(pheno, cols);
Net = BayesNet(nodes, '', [],[],[],true,[],data,cols,pheno,[],{},tree);
[acc, p, z] = PredictPheno(Net, false);
[testval, ~] = AUCWorker(acc,p,z,data(:,phencol));

tpval = adaptivepermtester(numsims, testval, 1, @adaptivepredworker, ...
    tree, nodes, data, cols, pheno);
