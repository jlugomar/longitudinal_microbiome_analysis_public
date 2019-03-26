function DemoAnalysis(data, cols, pheno, analysis_title, priorPrecision, bootslimit, nodelimit)
% DemoAnalysis(data, cols, analysis_title, priorPrecision, bootslimit, nodelimit)
%
% do testing and demo most functions on basically any dataset, that 
% includes a testing dataset.  Input is two datasets in data{1} and
% data{2}. Can optionally not include seconod dataset; in which case
% testing will be omitted.
%
% focuses testing on cross validation and bootstrapping of each of four
% different search algorithms: 1) K2, 2) Pheno-Centric, 3) Full/Exhaustive
% greedy search, 4) naive bayes.
%
%
% INPUT
% DATA: data sets : DATA{1} is training DATA{2} is testing.
% COLS: column names (variable names) for DATA, 
%   COLS{1}: matching DATA{1}
%   COLS{2}: matching DATA{2}
% PHENO: string of phenotype column in COLS
% ANALYSIS_TITLE: string providing some kind of mnemonic representation of
%   the analysis being done here
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% BOOTSLIMIT: number of boostrap realizations to do for each algorithm.
%   Optional.  Default = 25.
% NODELIMIT: maximum number of nodes to include in consensus networks 
%   computed from bootstrap networks. Optional.  Default = 25.
%
%
% copyright Michael McGeachie, 2013. MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    bootslimit = 25;
end
if (nargin < 7)
    nodelimit = 25;
end

if (iscell(data) && length(data) > 1)
    tdatas{1} = data{1};
    tdatas{2} = data{2};
    tcols{1} = cols{1};
    tcols{2} = cols{2};
    TESTING = true;
else
    TESTING = false;
    tdatas{1} = data;
    tcols{1} = cols;
end

verbose = false;

%% do K2
[K2MBnet, K2FullBN] = LearnStructure(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    [analysis_title,'K2'], verbose);
K2FullBN.WriteToGML([analysis_title,'K2_full']);
K2MBnet.WriteToGML([analysis_title,'K2_MB']);
auc_k2 = BNLearnAndTest(K2MBnet);
if (TESTING)
    fprintf(1,'Do testing of K2 algorithm on test data: \n');
    tauc_k2 = BNLearnAndTest(K2MBnet, tdatas{2},tcols{2});
    fprintf(1,'Testing performance of K2 algorithm AUC : %2.5f\n', 100 * tauc_k2);
end

%%
fprintf(1,'Do CV building of K2 algorithm\n');
% build different networks in CV (don't output them)
algorithm = 1;
folds = 5;
verbose = true;
[auc_k2cv, k2cv_numnodes, classacc, cvAdjMat] = CVModelLearnEx(tdatas{1}, tcols{1}, ...
    pheno, priorPrecision, folds, [analysis_title,'K2_CV'], verbose, false, algorithm);
verbose = false;
% retain the AUC in CV testing
fprintf('total K2 CV AUC : %2.5f\n',100 * auc_k2cv);

%% do CV, params only:
fprintf(1,'Do CV on paramters from K2 network\n');
folds = 5;
[auc_k2cvparams, ~] = CVParamTest(K2MBnet, folds, '', verbose, false);
fprintf('total K2 network CV (params) AUC : %2.5f\n',100 * auc_k2cvparams);

%% do bootstrapping to build a K2 network:
algorithm = 1;
fprintf(1,'Do Bootstrap building of K2 network\n');
BootsAdjMat = BootstrapLearn(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    bootslimit, algorithm, verbose);
if (TESTING)
    [k2_bsAUCs, k2_bsAUCs_test, k2_bsnumnodes] = BootsNetAddEdgesCVandTest(BootsAdjMat, ...
        tdatas{1}, tcols{1}, tdatas{2}, tcols{2}, priorPrecision, pheno, nodelimit);
else
    [k2_bsAUCs, k2_bsAUCs_test, k2_bsnumnodes] = BootsNetAddEdgesCV(BootsAdjMat, ...
        tdatas{1}, tcols{1},  priorPrecision, pheno, nodelimit);    
end


%% do phenocentric
BF_THRESH = 0;
[PCMBnet] = LearnPhenoCentric(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    BF_THRESH, verbose);
PCMBnet.title = [analysis_title,'PC_MB'];
PCMBnet.WriteToGML();
auc_pc = BNLearnAndTest(PCMBnet);
if (TESTING)
    fprintf(1,'Do testing of Pheno-Centric algorithm on test data: \n');
    tauc_pc = BNLearnAndTest(PCMBnet, tdatas{2},tcols{2});
    fprintf(1,'Testing performance of Pheno-Centric algorithm AUC : %2.5f\n', 100 * tauc_pc);
end

%%
fprintf(1,'Do CV building of PhenoCentric algorithm\n');
% build different networks in CV (don't output them)
algorithm = 2;
folds = 5;
verbose = true;
[auc_pccv, pccvnumnodes, classacc, cvAdjMat] = CVModelLearnEx(tdatas{1}, tcols{1}, ...
    pheno, priorPrecision, folds, [analysis_title,'PC_CV'], verbose, false, algorithm);
verbose = false;
% retain the AUC in CV testing
fprintf('total PhenoCentric CV AUC : %2.5f\n',100 * auc_pccv);

%% do CV, params only:
fprintf(1,'Do CV on paramters from PhenoCentric network\n');
folds = 5;
[auc_pccvparams, ~] = CVParamTest(PCMBnet, folds, '', verbose, false);
fprintf('total PhenoCentric network CV (params) AUC : %2.5f\n',100 * auc_pccvparams);


%% do bootstrapping to build a PhenoCentric network:
algorithm = 2;
fprintf(1,'Do Bootstrap building of PhenoCentric network\n');
BootsAdjMat = BootstrapLearn(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    bootslimit, algorithm, verbose);
if (TESTING)
    [pc_bsAUCs, pc_bsAUCs_test, pc_bsnumnodes] = BootsNetAddEdgesCVandTest(BootsAdjMat, ...
        tdatas{1}, tcols{1}, tdatas{2}, tcols{2}, priorPrecision, pheno, nodelimit);
else
    [pc_bsAUCs, pc_bsAUCs_test, pc_bsnumnodes] = BootsNetAddEdgesCV(BootsAdjMat, ...
        tdatas{1}, tcols{1}, priorPrecision, pheno, nodelimit);
end

%% do naive bayes (set maxParents = 1):
BF_THRESH = 0;
priorPrecision.maxParents = 1;
[NBnet] = LearnPhenoCentric(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    BF_THRESH, verbose);
NBnet.title = [analysis_title,'NaiveBayes'];
NBnet.WriteToGML();
auc_nb = BNLearnAndTest(NBnet);
if (TESTING)
    fprintf(1,'Do testing of Naive Bayes algorithm on test data: \n');
    tauc_nb = BNLearnAndTest(NBnet, tdatas{2},tcols{2});
    fprintf(1,'Testing performance of Naive Bayes algorithm AUC : %2.5f\n', 100 * tauc_nb);
end
priorPrecision.maxParents = 2;

%%
fprintf(1,'Do CV building of Naive Bayes algorithm\n');
% build different networks in CV (don't output them)
algorithm = 2;
folds = 5;
verbose = true;
priorPrecision.maxParents = 1;
[auc_nbcv, numnodes, classacc, cvAdjMat] = CVModelLearnEx(tdatas{1}, tcols{1}, ...
    pheno, priorPrecision, folds, [analysis_title,'NB_CV'], verbose, false, algorithm);
priorPrecision.maxParents = 2;
verbose = false;
% retain the AUC in CV testing
fprintf('total PhenoCentric CV AUC : %2.5f\n',100 * auc_nbcv);

%% do CV, params only:
fprintf(1,'Do CV on paramters from Naive Bayes network\n');
folds = 5;
priorPrecision.maxParents = 1;
[auc_nbcvparams, ~] = CVParamTest(NBnet, folds, '', verbose, false);
priorPrecision.maxParents = 2;
fprintf('total Naive Bayes network CV (params) AUC : %2.5f\n',100 * auc_nbcvparams);

%% do bootstrapping to build a Naive Bayes network:
algorithm = 2;
priorPrecision.maxParents = 1;
fprintf(1,'Do Bootstrap building of Naive Bayes network\n');
BootsAdjMat = BootstrapLearn(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    bootslimit, algorithm, verbose);
if (TESTING)
    [nb_bsAUCs, nb_bsAUCs_test, nb_bsnumnodes] = BootsNetAddEdgesCVandTest(BootsAdjMat, ...
        tdatas{1}, tcols{1}, tdatas{2}, tcols{2}, priorPrecision, pheno, nodelimit);
else
    [nb_bsAUCs, nb_bsAUCs_test, nb_bsnumnodes] = BootsNetAddEdgesCV(BootsAdjMat, ...
        tdatas{1}, tcols{1}, priorPrecision, pheno, nodelimit);    
end
priorPrecision.maxParents = 2;

%% do full netowrk search
BF_THRESH = 0;
[FullBNet] = FullBNLearn(tdatas{1}, tcols{1}, pheno, BF_THRESH, [analysis_title,'Full'], priorPrecision);
FullBN_MB = FullBNet.MakeIntoMB();
FullBNet.title = [analysis_title,'Full_full'];
FullBN_MB.title = [analysis_title,'Full_MB'];
FullBNet.WriteToGML();
FullBN_MB.WriteToGML();
auc_full = BNLearnAndTest(FullBN_MB);
if (TESTING)
    fprintf(1,'Do testing of Full-Exhaustive algorithm on test data: \n');
    tauc_full = BNLearnAndTest(FullBN_MB, tdatas{2},tcols{2});
    fprintf(1,'Testing performance of Full-Exhaustive algorithm AUC : %2.5f\n', 100 * tauc_full);
end

%%
fprintf(1,'Do CV building of Full-Exhaustive algorithm\n');
% build different networks in CV (don't output them)
algorithm = 3;
folds = 5;
verbose = true;
[auc_fullcv, numnodes, classacc, cvAdjMat] = CVModelLearnEx(tdatas{1}, tcols{1}, ...
    pheno, priorPrecision, folds, [analysis_title,'Full_CV'], verbose, false, algorithm);
verbose = false;
% retain the AUC in CV testing
fprintf('total PhenoCentric CV AUC : %2.5f\n',100 * auc_fullcv);

%% do CV, params only:
fprintf(1,'Do CV on paramters from Full network\n');
folds = 5;
[auc_fullcvparams, ~] = CVParamTest(FullBN_MB, folds, '', verbose, false);
fprintf('total Full network CV (params) AUC : %2.5f\n',100 * auc_fullcvparams);

%% do bootstrapping to build a Full-Exhaustive network:
algorithm = 3;
fprintf(1,'Do Bootstrap building of Full-Exhaustive network\n');
BootsAdjMat = BootstrapLearn(tdatas{1}, tcols{1}, pheno, priorPrecision, ...
    bootslimit, algorithm, verbose);
if (TESTING)
    [full_bsAUCs, full_bsAUCs_test, full_numnodes] = BootsNetAddEdgesCVandTest(BootsAdjMat, ...
        tdatas{1}, tcols{1}, tdatas{2}, tcols{2}, priorPrecision, pheno, nodelimit);
else
    [full_bsAUCs, full_bsAUCs_test, full_numnodes] = BootsNetAddEdgesCV(BootsAdjMat, ...
        tdatas{1}, tcols{1}, priorPrecision, pheno, nodelimit);    
end


%% make a nice graph of the bootstrap results

figure();
hold on;
m = min(length(k2_bsnumnodes),length(k2_bsAUCs));
plot(k2_bsnumnodes(2:m), k2_bsAUCs(2:m), '-b*');
m = min(length(k2_bsnumnodes),length(k2_bsAUCs_test));
plot(k2_bsnumnodes(2:m), k2_bsAUCs_test(2:m),':bo');
m = min(length(pc_bsnumnodes),length(pc_bsAUCs));
plot(pc_bsnumnodes(2:m), pc_bsAUCs(2:m), '-r*');
m = min(length(pc_bsnumnodes),length(pc_bsAUCs_test));
plot(pc_bsnumnodes(2:m), pc_bsAUCs_test(2:m),':ro');
m = min(length(nb_bsnumnodes),length(nb_bsAUCs));
plot(nb_bsnumnodes(2:m), nb_bsAUCs(2:m), '-k*');
m = min(length(nb_bsnumnodes),length(nb_bsAUCs_test));
plot(nb_bsnumnodes(2:m), nb_bsAUCs_test(2:m),':ko');
m = min(length(full_numnodes),length(full_bsAUCs));
plot(full_numnodes(2:m), full_bsAUCs(2:m), '-m*');
m = min(length(full_numnodes),length(full_bsAUCs_test));
plot(full_numnodes(2:m), full_bsAUCs_test(2:m),':mo');
hold off;
set(gca,'Title',text('String','AUC in CV and Test for Included Nodes from Bootsrap Networks'));
xlabel('Number of nodes in Markov Blanket');
ylabel('AUC %');
m = max([k2_bsnumnodes, pc_bsnumnodes, nb_bsnumnodes, full_numnodes]);
axis([1,m,.5,1]);
legend('K2 CV', 'K2 Test', 'PhenoCentric CV', 'PhenoCentric Test', ...
    'NaiveBayes CV', 'NaiveBayes Test', 'FullSearch CV','FullSearch Test'); 
%
% writeFig300dpi(1, 'BootstrapMethodsCompare_Geo1.tiff');
%




