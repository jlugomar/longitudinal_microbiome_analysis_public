% script: TestSuite
% testing file for hybrid bayesian networks
%
% Loads various types of networks, then learns and predicts on them.
% Compare output to the output provided in the file "testsuite_output.txt"
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

%close all; clear all;
fprintf(1,'\n\n======  Test Bayes Net Inference =====\n\n');
run('../../bnpathscript');

% timing
tic;

% load a network of Cowell's Figure 2 and figure 4.  
% It loads a network called 'cowell' and another called 'cfig4'
run('mnetscripts');

% common parameter values:
priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 3;

nodesketch = cowell;
[data,cols] = ReadHeaderDataFile('cowell.txt');
cowell_BN = BayesNet(nodesketch, 'cowell_fig2', [],[],[],[],[],data,cols,'A',priorPrecision);
cowell_fig2_auc = BNLearnAndTest(cowell_BN);
fprintf(1,'Cowell''s Figure 2 AUC: %2.1f\n',100 * cowell_fig2_auc);
%DoTestNetworkEval(nodesketch, 'cowell.txt', 'A', priorPrecision, 'Cowell''s Figure 2', true, false);

nodesketch = cfig4;
[data,cols] = ReadHeaderDataFile('cfig4.txt');
cowell_BN = BayesNet(nodesketch, 'cowell_fig4', [],[],[],[],[],data,cols,'A',priorPrecision);
cowell_fig4_auc = BNLearnAndTest(cowell_BN);
fprintf(1,'Cowell''s Figure 4 AUC: %2.1f\n',100 * cowell_fig4_auc);
%DoTestNetworkEval(nodesketch, 'cfig4.txt', 'A', priorPrecision, 'Cowell''s Figure 4', true, false);

% common parameter values:
priorPrecision.nu = 25;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 25;
priorPrecision.maxParents = 3;

[winedata,winecols] = ReadHeaderDataFile('winedata.txt');
% test BFFilterBNLearn with various algrothms:
fprintf(1,'\nLearning Wine example with K2 Search\n');
[auc,MBNet] = BFFilterBNLearn(winedata,winecols,'class',1,0,true,priorPrecision);
MBNet.title = 'GMLTest1-K2';
MBNet.WriteToGML();
fprintf(1,'-- Wine K2: %2.1f AUC on %d nodes\n', 100*auc, length(MBNet.mb));
fprintf(1,'Learning Wine example with Pheno-Centric Search\n');
[auc,MBNet] = BFFilterBNLearn(winedata,winecols,'class',2,0,true,priorPrecision);
MBNet.title = 'GMLTest2-PhenoCentric';
MBNet.WriteToGML();
fprintf(1,'-- Wine PhenoCentric: %2.1f AUC on %d nodes\n', 100*auc, length(MBNet.mb));
fprintf(1,'Learning Wine example with Exhaustive Search\n');
[auc,MBNet] = BFFilterBNLearn(winedata,winecols,'class',3,0,true,priorPrecision);
MBNet.title = 'GMLTest3-Exhaustive';
MBNet.WriteToGML();
fprintf(1,'-- Wine Exhaustive: %2.1f AUC on %d nodes\n', 100*auc, length(MBNet.mb));


% test Bootstrap learning using WINE:
fprintf(1,'\nDoing K2 Bootstrap Testing on WINE dataset\n');
BootsAdjMat = BootstrapLearn(winedata, winecols, 'class', priorPrecision, 5, 1);
% then make a bootstrap network, and test it:
[WBN, WineBootsNetMB] = NEdgesFromBoots(BootsAdjMat, winedata, winecols, 'class', 6);
WBN.WriteToTGF();
% measure AUC
WineBootsNetMB.priorPrecision = priorPrecision;
[auc, p, z] = BNLearnAndTest(WineBootsNetMB);
fprintf(1,'-- Wine K2 Boots: %2.1f AUC on %d nodes\n', 100*auc, length(WineBootsNetMB.mb));

BootsAdjMatPhenoCentric = BootstrapLearn(winedata, winecols, 'class', priorPrecision, 5, 2);
[WBN, WineBootsNetMB] = NEdgesFromBoots(BootsAdjMatPhenoCentric, winedata, winecols, 'class', 6);
WBN.WriteToTGF();
% measure AUC
WineBootsNetMB.priorPrecision = priorPrecision;
[auc, p, z] = BNLearnAndTest(WineBootsNetMB);
fprintf(1,'-- Wine PhenoCentric Boots: %2.1f AUC on %d nodes\n', 100*auc, length(WineBootsNetMB.mb));

BootsAdjMatFull = BootstrapLearn(winedata, winecols, 'class', priorPrecision, 5, 3);
[WBN, WineBootsNetMB] = NEdgesFromBoots(BootsAdjMatFull, winedata, winecols, 'class', 6);
WBN.WriteToTGF();
% measure AUC
WineBootsNetMB.priorPrecision = priorPrecision;
[auc, p, z] = BNLearnAndTest(WineBootsNetMB);
fprintf(1,'-- Wine Exhaustive Boots: %2.1f AUC on %d nodes\n', 100*auc, length(WineBootsNetMB.mb));



%% test learning a BN and it's prediction on WINE.
[data, cols] = ReadHeaderDataFile('winedata.txt');
MBNet = LearnStructure(data, cols, 'class', priorPrecision, 'winetest-net');
[auc, p, z] = BNLearnAndTest(MBNet);
%DoTestNetworkEval(MBNet, 'winedata.txt', 'class', priorPrecision, ...
%    'UCI Machine Learning Repository: WINE', false, false);

% test inference on a fully discrete bayes net
% also test the input/output to BayesWare Discoverer
[acc,p] = BWDHybridInference('chess-kr-vs-kp.bdn', 'chess-kr-vs-kp.txt', 'class');

% now learn the chess example using LearnStructure:
[data, cols] = ReadHeaderDataFile('chess-kr-vs-kp.txt');
MBNet = LearnStructure(data,cols,'class',priorPrecision, 'chess-net');
[auc, p, z] = BNLearnAndTest(MBNet);

%% ModelLearnAndTest
% should test this with a dataset with many discrete variables:
[data, cols] = ReadHeaderDataFile('winedata.txt');
pheno = 'class';
verbose = true;
rvec = rand(1,length(data(:,1)));
traindata = data(rvec > 0.5,:);
testdata = data(rvec < 0.5,:);
% common parameter values:
priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 2;
fprintf(1,'-- Running ModelLearnAndTest on Wine Data with Random Testing/Training\n');
algorithm = 1;
[auc, numnodes, testauc, model] = ModelLearnAndTest(traindata, cols, testdata, cols, ...
    pheno, priorPrecision, 'chess-ModelLearnAndTest-K2', verbose, {}, algorithm);

algorithm = 2;
[auc, numnodes, testauc, model] = ModelLearnAndTest(traindata, cols, testdata, cols, ...
    pheno, priorPrecision, 'chess-ModelLearnAndTest-PhenoCentric', verbose, {}, algorithm);

algorithm = 3;
[auc, numnodes, testauc, model] = ModelLearnAndTest(traindata, cols, testdata, cols, ...
    pheno, priorPrecision, 'chess-ModelLearnAndTest-Full', verbose, {}, algorithm);



%% do DemoAnalysis() with and without a testing set.
demodata = {traindata,testdata};
democols = {cols,cols};
DemoAnalysis(demodata, democols, pheno, 'chess-DemoAnalysis-Test', priorPrecision);



t = toc;

fprintf(1,' *** Total time Elapsed: %d ***\n',t);
