% function script:
% load and analyze human_cachexia dataset, publicly available metabolomic
% dataset from MetaboAnalyst(2.0).
run '../../bnpathscript';
% load dataset
[numdata, cols, strcols, numcolindex, strcolindex] = RCSVLoad('human_cachexia.csv', false);
pheno = 'Muscle loss';

% replace string pheno with binary pheno
cases = strcmp('cachexic', strcols{2});
phncol = zeros(size(strcols{2}));
phncol(cases) = 1;
data = [phncol, numdata];
cols = cols(2:end);

%% normalize
% first log xform all metabolites in data2:
FIRSTMET = 2; % looks like first metabolite is in the 2nd column
for i = FIRSTMET:size(data,2)
    data(:,i) = log2(data(:,i));
end

disc = IsDiscrete(data);
for i = 1:length(disc)
    if (~disc(i))
        % normalize:
        md = mean(data(:,i));
        stdd = std(data(:,i));
        data(:,i) = (data(:,i) - md) ./ stdd;
    end
end

%% do a CV network parameter exploration
% probably skip this due to low sample size; variance in AUCs from each
% hold out set is too large for reasonable results.
%analysis_title = 'Human Cachexia Demo';
% time consuming step:
%[bestParams,aucs] = CVExploreParams(data, cols, pheno, analysis_title);

%% learn a network
analysis_title = 'Human Cachexia Demo';
fprintf(1,'\n\n======  start new analysis: %s  =====\n\n',analysis_title);
tic;

% common parameter values:
priorPrecision.nu = 10;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 2;

verbose = true;
doAUC = true;
fprintf(1,'Learning Network Structure for %s\n', analysis_title);
MBNet = LearnStructure(data, cols, pheno, priorPrecision, [analysis_title,'-net'], verbose);
fprintf(1,'Learning Network Parameters\n');
MBNet = LearnParamsBN(MBNet);
fprintf(1,'Predicting on Training Data\n');
[acc, p, z] = PredictPheno(MBNet, verbose);

% output AUC or accuracy measure:
convexhullAUC = AUCWorker(acc,p,z,data(:,1),doAUC);


%% learn a bootstrapped network
%nboots = 250;
nboots = 25;
cadjmat = BootstrapLearn(data, cols, pheno, priorPrecision, nboots);

% compute networks at various edge-frequency thresholds.
%thresholds = [.20, .25, .30, .50];
%[bootAUCs, bootpvps] = ThreshBootNet(cadjmat, data, cols, pheno, priorPrecision, verbose, thresholds);

%% view bootstrap edges in order of prediction accuracy impact:
[cvAUCs, ~, ~] = BootsNetAddEdgesCV(cadjmat, ...
    data, cols, priorPrecision, pheno);
figure();
plot(cvAUCs);
title('AUC by number of nodes included in BN.');
xlabel('Number of Nodes, ordered by most common in Bootstrap networks');
ylabel('AUC in 5-fold Cross Validation');

%%
% make 4 models which make different tradeoffs between model complexity and 
% model performance. 
numnodes = [8,4,3,2];
for i = 1:length(numnodes)
    [BN, TMBNet{i}] = NEdgesFromBoots(cadjmat, data, cols, pheno, numnodes(i));
    outputBayesNetGraph(TMBNet{i}.adjmat,TMBNet{i}.cols, ['CachexiaNet', num2str(numnodes(i)),'Nodes.tgf']);
    
    % measure AUC
    TMBNet{i}.nodes = BNfromAdjMat(TMBNet{i}.adjmat, TMBNet{i}.disc, TMBNet{i}.cols);
    TMBNet{i}.priorPrecision = priorPrecision;
    TMBNet{i} = LearnParamsBN(TMBNet{i});
    [accs, ps, zs{i}] = PredictPheno(TMBNet{i}, false);
    bootnodes_auc(i) = AUCWorker(accs, ps, zs{i}, TMBNet{i}.data(:,1), true, false, true, true);
end

% get pvalues for difference between these four networks
pvp = -1*ones(length(numnodes));
for i = 1:length(numnodes)-1
    for j = i+1: length(numnodes)
        pvp(i,j) = CompTwoAUC(data(:,1), zs{i}, data(:,1), zs{j}, true);
    end
end

