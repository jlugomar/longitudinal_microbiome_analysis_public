function tpval = adaptivecvworker(adjmat, data, cols, pheno, prior, discrete, folds)
%
% DEPRICATED
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

[ncases,~] = size(data);
phencol = strmatch(pheno, cols, 'exact');

% just permute the labels (or resample them):
r = randperm(ncases);
simdata = data;
simdata(:,phencol) = simdata(r,phencol);

% split data into N-fold CV sets:
r = ceil(rand(1,ncases) * folds);

num = zeros(folds, 1);
acc = num;
cvPClass = [];
cvPredZs = [];
cvTrueClass = [];

% check values of discrete vars::
discvals = cell(size(discrete));
for i = 1:length(cols)
    if (discrete(i))
        discvals{i} = num2cell(unique(data(:,i),'legacy'));
    else 
        discvals{i} = {};
    end
end

% do a CV-splits
for i = 1:folds
    cvdata = simdata(r ~= i,:);
    cvtest = simdata(r == i,:);
    if (isempty(cvtest))
        continue;
    end
    num(i) = length(cvtest);

    MBNet = BNfromAdjMat(adjmat, discrete, cols);
    [tree, nodes] = LearnParams(MBNet , '', cvdata, cols, prior, discvals, discrete);
    [acc(i), p, z] = PredictPheno(tree, nodes, '', pheno, cvtest, cols, false);

    cvPClass = [cvPClass; p];
    cvPredZs = [cvPredZs; z];
    cvTrueClass = [cvTrueClass; cvtest(:,phencol)];
end

acc_all = acc * num' ./ sum(num);

[tpval, ~] = AUCWorker(acc_all,cvPClass,cvPredZs,cvTrueClass);
