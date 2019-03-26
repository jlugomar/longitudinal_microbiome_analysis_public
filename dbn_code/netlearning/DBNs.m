% dynamic bayesian networks
% plan:
% 1) combine timeseries data into a larger column matrix with each time point
% matrix below the prior time point matrix.
% 2) learn network using StateTrackerSearch with dynamic bayes net option
% enabled
% ... allow cycles, allow self-loops
% 3) unroll BN tinto a 2TBN : a 2-timepoint BN; all edges are from time
% point one to time point two.
% 3a) unroll dataset from timeseries matrix
% 4) use normal techniques to predict with unrolled 2TBN
%
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

% do a simple test of a DBN:

%% 1) combine Geo dataset :
% /data/testdata/Geo/
[data,subjids,cols] = ProcessGeoFor2TDBN();

%% compute BayesFactors and take a subset of the data:
pheno = 'Exacerbation';
disc = IsDiscrete(data);
verbose = false;
priorPrecision.nu = 10;
priorPrecision.alpha = 10;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = 3;

bfs = BayesFactorScore(data, cols, pheno, priorPrecision, disc, verbose);
keep = bfs > 10;
% also keep the first 9 columns : these are some clinical ones and the
% phenotypes.
keep(1:9) = true;
data = data(:,keep);
cols = cols(keep);

% drop 'Donor' column which is just the subject id:
dropcols = {'Donor'};
drops = false(size(cols));
for i = 1:length(dropcols)
    m = strcmp(dropcols(i),cols);
    drops = drops | m;
end
data = data(:,~drops);
cols = cols(~drops);


%% take data and split into two (overlapping) time points:
[d0,dn] = MakeTSBNData(data, subjids);



%% make main call to learn the DBN:
BFTHRESH = 0;
outfilename = 'DBN_GEO_test1';
disc = IsDiscrete(data);
searchParameter.t0cont = d0(:,~disc)';
searchParameter.tncont = dn(:,~disc)';
searchParameter.t0disc = d0(:,disc)';
searchParameter.tndisc = dn(:,disc)';
searchParameter.DBN = true;
verbose = true;
[BNet, outstats] = FullBNLearn(data, cols, pheno, BFTHRESH, outfilename, ...
    priorPrecision, disc, verbose, searchParameter);

%% reduce to DBN network, then print out:
DBN = BNet.ConvertToDBN();
DBN.WriteToGML('DBN_testv1_DBN');
DBN.WriteToTGF('DBN_testv1_DBN');









