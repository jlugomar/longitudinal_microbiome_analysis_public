%% test with GSE19301_diffexpr.txt
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

fprintf(1,'Loading and processing GSE190301 differential expression data\n');
[data, cols, strdata, dcolinds, scolinds] = RCSVLoad('GSE19301_diffexpr_geneids.txt',false,'\t');
% make exacerbation under greater steroid exposure the phenotype:
pheno = 'maximum steroid exposure  4=systemic  3=inhaled  2=intranasal  1';

dropcols = {'Exacerbation','Donor','leukotriene receptor antagonist', ...
    'baseline severity nih guideline'};

drops = zeros(size(cols));
for i = 1:length(dropcols)
    drops = drops | strcmpi(dropcols(i),cols);
end

% drop the columns we're not interested in
cols = cols(~drops);
data = data(:,~drops);

% move phenotype to column #1
phenocol = strcmpi(pheno,cols);
data = [data(:,phenocol),data(:,~phenocol)];
% and rename phenotype:
pheno = 'Steroid Exacerbator';
cols = {pheno, cols{~phenocol}};

% discretize phenotype: (indicator of someone that has an exacerbation
% while on more steroids that usual)
data(:,1) = data(:,1) > 0;
% set discrete indicator
disc = [true, false(1,length(cols)-1)];
RACE = 'race';
racecol = find(strcmpi(RACE,cols));
disc(racecol) = true;

% normalize all other datacolumns:
fprintf(1,'Normalizing continuous data columns\n');
for i=1:length(disc)
    if (~disc(i))
        data(:,i) = (data(:,i) - mean(data(:,i))) / std(data(:,i));
    end
end



%% now test some bayesian stuff:
priorPrecision.nu = 50;
priorPrecision.alpha = 50;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = 2;

% test bayes factors: BFScore
fprintf(1,'Computing Bayes Factors for all %d variables\n',length(cols));
BFs = BayesFactorScore(data,cols,pheno,priorPrecision);
fprintf(1,'Found %d variables with positive Bayes Factor\n',sum(BFs > 0));

%% do Naive Bayes Net w/o splitting data:
BFTHRESH = 0;
tic;
fprintf(1,'Building Naive Bayes model over %d variables\n',length(cols));
verbose = true;
NBBN = NBBayesNet(data, cols, pheno, BFTHRESH, 'Geo_full_NB', priorPrecision, disc, verbose);
t = toc;
fprintf(1,'Done Building Naive Bayes BN in %5.2f seconds\n',t);
NBBN.WriteToGML('GEO_full_NaiveBayes');

% test NB looking at full data to check fit for prediction:
fprintf(1,'Learning network distribution parameters\n');
tic;
NBBN = NBBN.LearnParams();
fprintf(1,'Predicting on with NB Model on all data\n');
% slow:
[acc, p, z] = NBBN.Predict(verbose);
[auc, classacc] = AUCWorker(acc,p,z,NBBN.GetPhenoCol());
% gets 83% accuracy and 89% AUC.
t = toc;
fprintf(1,'Done predicting NB Model on all data in %5.2f seconds\n',t);


%%
% normally, a random split is advised.  However, for purposes of
% reproducibility, just take the first half:
%rvec = rand(1,length(data(:,1)));
splitvec = zeros(length(data(:,1)),1);
splitvec(length(splitvec)/2+1:end) = 1;
traindata = data(splitvec > 0.5,:);
testdata = data(splitvec < 0.5,:);

% try Naive Bayes:
BFTHRESH = 0;
tic;
fprintf(1,'Building Naive Bayes model over %d variables\n',length(cols));
verbose = true;
NBBN = NBBayesNet(traindata, cols, pheno, BFTHRESH, 'NaiveBayesNetwork', priorPrecision, disc, verbose);
t = toc;
fprintf(1,'Done Building Naive Bayes BN in %5.2f seconds\n',t);
NBBN.WriteToGML('GEO_test_NaiveBayes');

%%
% test NB
fprintf(1,'Learning network distribution parameters\n');
NBBN = NBBN.LearnParams();
%
fprintf(1,'Predicting on Test Data\n');
NBBNtest = NBBN.ReplaceData(testdata,cols);
% very slow:
[acc, p, z] = NBBNtest.Predict(verbose);
% accuracy = 57% which should be significant
%
[auc, classacc] = AUCWorker(acc,p,z,NBBNtest.GetPhenoCol());

% check this AUC with a p-value:
numsims = 1000;
[permpval] = PermTestRandomAUC(auc, testdata, cols, pheno, numsims);


%% filter for everything else:
% finds ~250 variables. 
keep = BFs > 3;
% keep race in the mix:
racecol = find(strcmpi(RACE,cols));
keep(racecol) = true;
keep(1) = true;

cols = cols(keep);
data = data(:,keep);

% normally, a random split is advised.  However, for purposes of
% reproducibility, just take the first half:
%rvec = rand(1,length(data(:,1)));
splitvec = zeros(length(data(:,1)),1);
splitvec(length(splitvec)/2+1:end) = 1;
traindata = data(splitvec > 0.5,:);
testdata = data(splitvec < 0.5,:);


%% Try building a Pheno-Centric BN:
% phenocentric search

BFTHRESH = 0;
tic;
fprintf(1,'Building PhenoCentric BN over %d variables\n',length(cols));
PhenoCentricBN = LearnPhenoCentric(data, cols, pheno, priorPrecision, BFTHRESH, true, disc);
t = toc;
fprintf(1,'Done Building PhenoCentric BN in %5.2f seconds\n',t);
PhenoCentricBN.WriteToGML('GEO_test_phenocentric_full');


%% K2 search
tic;
fprintf(1,'Building K2 BN over %d variables\n',length(cols));
[K2MBnet, K2FullBN] = LearnStructure(data, cols, pheno, priorPrecision, 'K2NetTest', true);
t = toc;
fprintf(1,'Done Building K2 BN in %5.2f seconds\n',t);
K2FullBN.WriteToGML('GEO_test_K2_full');



%% exhaustive search

BFTHRESH = 0;
fprintf(1,'Building Exhaustive BN over %d variables\n',length(cols));
tic;
ExhaustiveBN = FullBNLearn(data, cols, pheno, BFTHRESH, 'ExhaustiveBN', priorPrecision, disc);
t = toc;
fprintf(1,'Done Building Exhaustive BN in %5.2f seconds\n',t);
ExhaustiveBN.WriteToGML('GEO_test_Exhaustive_full');


%% ModelLearnAndTest
% should test this with many discrete vars (in TestSuite.m):
algorithm = 1;
verbose = true;
[auc, numnodes, testauc, model, classacc, testclassacc] = ...
    ModelLearnAndTest(traindata, cols, testdata, cols, pheno, priorPrecision, ...
    'MLearnAndTest1', verbose, {}, algorithm);


%% BFFilterBNLearn()
BFTHRESH = 6;
algorithm = 1;
[auc,MBNet] = BFFilterBNLearn(data, cols, pheno, algorithm, BFTHRESH, verbose, priorPrecision);
algorithm = 2;
[auc,MBNet] = BFFilterBNLearn(data, cols, pheno, algorithm, BFTHRESH, verbose, priorPrecision);
algorithm = 3;
[auc,MBNet] = BFFilterBNLearn(data, cols, pheno, algorithm, BFTHRESH, verbose, priorPrecision);

%% Run additional Demonstration functions in DemoAnalysis()
% Test CV
% Test Bootstrapping
tdatas = {traindata,testdata};
tcols = {cols,cols};
DemoAnalysis(tdatas, tcols, pheno, 'GeoTestDemo', priorPrecision, 5, 20);

%% test ModelLearnAndTest
verbose = true;
algorithm = 1;
[auc, numnodes, testauc, BNModel, classacc, testclassacc] = ...
    ModelLearnAndTest(traindata, cols, testdata, cols, pheno, priorPrecision, ...
    'MLAT_GeoTest', verbose, {}, algorithm);
algorithm = 2;
[auc, numnodes, testauc, BNModel, classacc, testclassacc] = ...
    ModelLearnAndTest(traindata, cols, testdata, cols, pheno, priorPrecision, ...
    'MLAT_GeoTest', verbose, {}, algorithm);
algorithm = 3;
[auc, numnodes, testauc, BNModel, classacc, testclassacc] = ...
    ModelLearnAndTest(traindata, cols, testdata, cols, pheno, priorPrecision, ...
    'MLAT_GeoTest', verbose, {}, algorithm);

%% LogFitAndTest
verbose = true;
[auc, good, numvars, testauc, model, classacc, testclassacc] = ...
    LogFitAndTest(traindata, cols, testdata, cols , pheno, verbose);

%% RandomNetworkPVal()
fprintf(1,' --- Testing capacity for overfitting using RandomNetworkPVal ---\n');
algorithm = 1;
numsims = 10;
[tpval] = RandomNetworkPVal(priorPrecision, data, cols, pheno, numsims, algorithm);
fprintf(1,'Empirical p-value from label permutation testing: %f\n',tpval);


%% ThreshBootNetAndTest()
% do some phenocentric boots nets:

algorithm = 2; % pheno-centric search
verbose = true;
nboots = 25;
bootsadjmat = BootstrapLearn(traindata, cols, pheno, priorPrecision, nboots, ...
    algorithm, verbose);

thresh = 0.5;
convexhullAUC = ThreshBootNetAndTest(bootsadjmat, traindata, cols, testdata, ...
    cols, pheno, priorPrecision, verbose, thresh);

fprintf(1,' --- test a 10-node Bootstrap Network \n');
nnodes = 10;
[BN, MBnet] = NEdgesFromBoots(bootsadjmat, traindata, cols, pheno, nnodes);
MBnet.priorPrecision = priorPrecision;
auc = BNLearnAndTest(MBnet, traindata, cols);
auc = BNLearnAndTest(MBnet, testdata, cols);
fprintf(1,'\t 10-node Bootstrap Network AUC : %2.1f \n',auc * 100);
MBnet.WriteToGML('GEO_test_10nodePhenoCentricBoots');

fprintf(1,' --- test a 25-node Bootstrap Network \n');
nnodes = 25;
[BN, MBnet] = NEdgesFromBoots(bootsadjmat, traindata, cols, pheno, nnodes);
MBnet.priorPrecision = priorPrecision;
auc = BNLearnAndTest(MBnet, traindata, cols);
auc = BNLearnAndTest(MBnet, testdata, cols);
fprintf(1,'\t 25-node Bootstrap Network AUC : %2.1f \n',auc * 100);
MBnet.WriteToGML('GEO_test_25nodePhenoCentricBoots');

fprintf(1,' --- test a 100-node Bootstrap Network \n');
nnodes = 100;
[BN, MBnet] = NEdgesFromBoots(bootsadjmat, traindata, cols, pheno, nnodes);
MBnet.priorPrecision = priorPrecision;
auc = BNLearnAndTest(MBnet, traindata, cols);
auc = BNLearnAndTest(MBnet, testdata, cols);
fprintf(1,'\t 100-node Bootstrap Network AUC : %2.1f \n',auc * 100);
MBnet.WriteToGML('GEO_test_100nodePhenoCentricBoots');
