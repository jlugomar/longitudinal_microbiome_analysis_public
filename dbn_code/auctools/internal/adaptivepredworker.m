function tpval = adaptivepredworker(tree, nodes, data, cols, pheno)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

[ncases,~] = size(data);

% just permute the labels (or resample them):
r = randperm(ncases);
simdata = data;
simdata(:,phencol) = simdata(r,phencol);
Net = BayesNet(nodes, '', [],[],[],true,[],simdata,cols,pheno,[],{},tree);
[acc, p, z] = PredictPheno(Net, false);
[tpval, ~] = AUCWorker(acc,p,z,Net.GetPhenoCol());

