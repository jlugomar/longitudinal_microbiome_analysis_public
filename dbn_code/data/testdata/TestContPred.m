% script: TestContPred
% testing file for continuous node prediction in hybrid bayesian networks
%
%
% Copyright Michael McGeachie, 2014.  MIT license. See cgbayesnets_license.txt.

%close all; clear all;
fprintf(1,'\n\n======  Test Bayes Net Inference =====\n\n');
run('../../bnpathscript');

% timing
tic;
% common parameter values:
priorPrecision.nu = 25;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 25;
priorPrecision.maxParents = 3;
BF_THRESH = 0;
searchParameter.annealing = false;

[winedata,winecols] = ReadHeaderDataFile('winedata_discex.txt');
outfilename = 'wine_FullSearch';
disc = IsDiscrete(winedata);
verbose = true;
% pick random cont node:
%pheno = 'alcalinity';
pheno = 'flavanoids';


%%
algorithms = {'Sim Annealing', 'Full Search', 'K2', 'Pheno-Centric'};

for i = 1:4
    fprintf('\nLearn %s with %s\n',pheno, algorithms{i});
    if (i == 1)
        searchParameter.annealing = true;
        [BNet, outstats] = FullBNLearn(winedata, winecols, pheno, BF_THRESH, outfilename, ...
            priorPrecision, disc, verbose);
    elseif (i == 2)
        searchParameter.annealing = false;
        [BNet, outstats] = FullBNLearn(winedata, winecols, pheno, BF_THRESH, outfilename, ...
            priorPrecision, disc, verbose,searchParameter);
    elseif (i == 3)
        [~,BNet] = LearnStructure(winedata, winecols, pheno, priorPrecision, 'K2_mbnet', verbose);
    elseif (i == 4)
        BNet = LearnPhenoCentric(winedata, winecols, pheno, priorPrecision, BF_THRESH, verbose);
    end
            % test using whole bayes net, instead of just markov blanket -
    % learn params
    fprintf('\nLearn parameters for %s network\n',algorithms{i});
    BNet = BNet.LearnParams();
    % predict the continuous phenotype:
    fprintf('\nPredict continuous phenotype with %s network\n',algorithms{i});
    [acc, p, z] = PredictPhenoCont(BNet, verbose);
    
    % if var is continuous, report error in avg z-score:
    truedata = BNet.GetPhenoCol(pheno);
    sd = std(truedata);
    zmiss = abs(p - truedata)./sd;
    fprintf('\tmean z-score error: %2.2f\n',mean(zmiss));
end






%% Partial Prediction:
% enter evidence on many of the vars; predict the other vars, or output
% their probabilities.
[BNet, outstats] = FullBNLearn(winedata, winecols, pheno, BF_THRESH, outfilename, ...
    priorPrecision, disc, verbose,searchParameter);
BNet = BNet.LearnParams();

basevars = 1:length(BNet.cols);
for i = 2:length(basevars)
    evinds = true(size(basevars));
    evinds(i) = false;
    evvars = basevars(evinds);
    predvars = [i];
    % includes continuous prediction
    fprintf(1,'Predicting variable %s in Wine network\n',BNet.cols{predvars});
    [preds, z, lep] = PredictPartial(BNet, evvars, predvars);
    
    % measure accuracy here:
    % if var is discrete, use AUC function
    if (disc(i))
        auc = AUCWorker(-1, preds, z, BNet.data(:,i),true,false);
        fprintf('\tAUC = %2.2f\n',auc);
    else
        % if var is continuous, report error in avg z-score:
        truedata = BNet.data(:,i);
        sd = std(truedata);
        zmiss = abs(preds - truedata)./sd;
        fprintf('\tmean z-score error: %2.2f\n',mean(zmiss));
    end
end

t = toc;

fprintf(1,' *** Total time Elapsed: %d ***\n',t);
