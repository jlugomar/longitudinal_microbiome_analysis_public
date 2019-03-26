function [t, acc] = OtestSuite()

% Parallel testing on Orchestra cluster
%
% testing file for hybrid bayesian networks
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

fprintf(1,'\n\n======  Test Bayes Net Inference =====\n\n');
addpath('../../inference/');
addpath('../../netlearning/');
addpath('../../auctools/');
addpath('../../utils/');

tic;

% load a network of Cowell's Figure 2 and figure 4.  
% It loads a network called 'cowell' and another called 'cfig4'
run('mnetscripts');

priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;

nodesketch = cowell;
DoTestNetworkEval(nodesketch, 'cowell.txt', 'A', priorPrecision, 'Cowell''s Figure 2', true, false);

nodesketch = cfig4;
DoTestNetworkEval(nodesketch, 'cfig4.txt', 'A', priorPrecision, 'Cowell''s Figure 4', true, false);

% test learning a BN and it's prediction on WINE.
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
%DoTestNetworkEval(MBNet, 'chess-kr-vs-kp.txt', 'class', priorPrecision, ...
%    'chess: Rook vs Pawn', true, false);

t = toc;

fprintf(1,' *** Total time Elapsed: %d ***',t);
