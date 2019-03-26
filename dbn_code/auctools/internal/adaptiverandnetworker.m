function tpval = adaptiverandnetworker(priorPrecision, data, cols, pheno, ALG, disc)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

[ncases,~] = size(data);
phencol = strmatch(pheno, cols, 'exact');

% just permute the labels (or resample them):
r = randperm(ncases);
simdata = data;
simdata(:,phencol) = simdata(r,phencol);

verbose = false;
BFTHRESH = 0;

if (ALG == 3)
    BN = FullBNLearn(simdata, cols, pheno, BFTHRESH, '', priorPrecision, disc);
    MBNet = BN.MakeIntoMB();
elseif (ALG == 2)
    MBNet = LearnPhenoCentric(simdata, cols, pheno, priorPrecision, BFTHRESH, verbose, disc);
else
    MBNet = LearnStructure(simdata, cols, pheno, priorPrecision, '', verbose, disc);
end

if (~isempty(MBNet) && length(MBNet.mb) > 1)
    [tpval] = BNLearnAndTest(MBNet, simdata, cols);
else
    tpval = 0.5;
end
