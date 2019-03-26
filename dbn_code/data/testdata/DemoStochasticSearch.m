
%% Test and Demo simulated annealing search procedure

%% first bit is just to load the Cachexia demo data:
% use the human_cachexia dataset, publicly available metabolomic
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

% normalize
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

%%  Actually do the stochastic searches:
analysis_title = 'Stochastic Search Procedures';
fprintf(1,'\n\n======  start new analysis: %s  =====\n\n',analysis_title);
tic;

% common parameter values:
priorPrecision.nu = 10;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 6;

verbose = true;
doAUC = true;
BF_THRESH = 0;
searchParameter.backtracking = true;
searchParameter.annealing = false;
fprintf(1,'Learning Network Structure w/ Exhaustive Hillclimber w/Backtracking for %s\n', analysis_title);
[ExBTR_Net, ex_stats1] = FullBNLearn(data, cols, pheno, BF_THRESH, [analysis_title,'-net'], ...
    priorPrecision, disc, verbose, searchParameter);
exbtr_lldiffs = ex_stats1.lldiffs;
exbtr_edges = ex_stats1.numedges;
exbtr_evals = ex_stats1.numevals;

searchParameter.backtracking = false;
searchParameter.annealing = false;
fprintf(1,'Learning Network Structure w/ Exhaustive Hillclimber NO Backtracking for %s\n', analysis_title);
[ExNo_Net, ex_stats2] = FullBNLearn(data, cols, pheno, BF_THRESH, [analysis_title,'-net'], ...
    priorPrecision, disc, verbose, searchParameter);
exno_lldiffs = ex_stats2.lldiffs;
exno_edges = ex_stats2.numedges;
exno_evals = ex_stats2.numevals;

%%
numruns = 10;
figure();
hold on;
c = colormap(jet);
[s1,~] = size(c);
ccstepsize = floor(s1 / (numruns + 4));
h = zeros (1,7);
h(1) = plot(exbtr_lldiffs, 'Color',c(ccstepsize * 1,:));
h(2) = plot(exno_lldiffs, 'Color',c(ccstepsize * 2,:));
plot(10 * exbtr_edges,':','Color',c(ccstepsize * 1,:));
plot(10 * exno_edges,':','Color',c(ccstepsize * 2,:));
hold off;

%% compare to phenocentric and K2 searches:

fprintf(1,'Learning Network Structure w/ Pheno Centric Search for %s\n', analysis_title);
[PC_MBNet, PC_FullNet, pc_outstats]=LearnPhenoCentric(data, cols, pheno, priorPrecision, BF_THRESH, verbose);

%% K2:

fprintf(1,'Learning Network Structure w/ K2 Search for %s\n', analysis_title);
[K2_MBnet, K2_FullBN, k2_outstats] = LearnStructure(data, cols, pheno, priorPrecision, [], verbose);

%% full and sim anneal:

BF_THRESH = 0;
priorPrecision.maxParents = 6;
searchParameter.backtracking = true;
searchParameter.annealing = false;
combo_lldiffs = cell(1,numruns);
combo_numedges = cell(1,numruns);
combo_temps = cell(1,numruns);
combo_numevals = cell(1,numruns);
for i = 1:numruns
    fprintf(1,'Learning Network Structure w/ Exhaustive Hillclimber followed by Sim Annealing for %s\n', analysis_title);
    [Combo_Net, combo_stats1] = FullAndSimAnneal(data, cols, pheno, BF_THRESH, [analysis_title,'-net'], ...
        priorPrecision, disc, verbose, searchParameter);
    combo_lldiffs{i} = combo_stats1{i}.lldiffs;
    combo_numedges{i} = combo_stats1{i}.numedges;
    combo_temps{i} = combo_stats1{i}.temp; 
    combo_numevals{i} = combo_stats1{i}.numevals;
end

%% do the simulated annealing stocastic search:
priorPrecision.maxParents = 6;
annealing = true;
backtracking = true;

SA_Net = cell(1,numruns);
sa_stats = cell(1,numruns);
lldiffs = cell(1,numruns);
numedges = cell(1,numruns);
temps = cell(1,numruns);
evals = cell(1,numruns);
searchParameter.backtracking = true;
searchParameter.annealing = true;
searchParameter.SA_Temp_Mult = 1;
for i = 1:numruns
    fprintf(1,'Learning Network Structure w/ Simulated Annealing for %s\n', analysis_title);
    [SA_Net{i}, sa_stats{i}] = FullBNLearn(data, cols, ...
        pheno, BF_THRESH, [analysis_title,'-net'], ...
        priorPrecision, disc, verbose, searchParameter);
    lldiffs{i} = sa_stats{i}.lldiffs;
    numedges{i} = sa_stats{i}.numedges;
    temps{i} = sa_stats{i}.temp; 
    evals{i} = sa_stats{i}.numevals;
end

%% plot the trajectories:
%figure();
hold on;
for i = 1:numruns
    h(3) = plot(lldiffs{i}, 'Color',c(ccstepsize * (i+2),:));
    h(4) = plot(10 * numedges{i},':','Color',c(ccstepsize * (i+2),:));
    h(5) = plot(100 * temps{i},'--','Color',c(ccstepsize * (i+2),:));
end
hold off;

title('Netowrk Search Strategies: Cachexia');
ylabel('Log Likelihood Improvement, Edges, & Temperature');
xlabel('Search Space (steps taken)');
%%

hold on;
h(6) = plot(k2_outstats.lldiffs, 'Color',c(ccstepsize * (numruns+3),:));
h(7) = plot(pc_outstats.lldiffs, 'Color',c(ccstepsize * (numruns+4),:));
plot(10 * k2_outstats.edges,':','Color',c(ccstepsize * (numruns+3),:));
plot(10 * pc_outstats.edges,':','Color',c(ccstepsize * (numruns+4),:));
hold off;
legend(h, 'Hill Climber', 'Backtracking Hill Cimber', 'Simulated Annealing', ...
    '10 * edges', 'SA Temperature', 'K2', 'Pheno-Centric');

%% display as a function of number of evals
numruns = 10;
figure();
hold on;
c = colormap(jet);
[s1,~] = size(c);
ccstepsize = floor(s1 / (numruns + 4));
h = zeros (1,5);
h(1) = plot(exbtr_evals, exbtr_lldiffs, 'Color',c(ccstepsize * 1,:));
h(2) = plot(exno_evals, exno_lldiffs, 'Color',c(ccstepsize * 2,:));
plot(exbtr_evals, 10 * exbtr_edges,':','Color',c(ccstepsize * 1,:));
plot(exno_evals, 10 * exno_edges,':','Color',c(ccstepsize * 2,:));
hold off;

hold on;
for i = 1:numruns
    h(3) = plot(evals{i}, lldiffs{i}, 'Color',c(ccstepsize * (i+2),:));
    h(4) = plot(evals{i}, 10 * numedges{i},':','Color',c(ccstepsize * (i+2),:));
    h(5) = plot(evals{i}, 100 * temps{i},'--','Color',c(ccstepsize * (i+2),:));
end
hold off;


hold on;
h(6) = plot(k2_outstats.evals, k2_outstats.lldiffs, 'Color',c(ccstepsize * (numruns+3),:));
h(7) = plot(pc_outstats.evals, pc_outstats.lldiffs, 'Color',c(ccstepsize * (numruns+4),:));
plot(k2_outstats.evals, 10 * k2_outstats.edges,':','Color',c(ccstepsize * (numruns+3),:));
plot(pc_outstats.evals, 10 * pc_outstats.edges,':','Color',c(ccstepsize * (numruns+4),:));
hold off;

legend(h, 'Hill Climber', 'Backtracking Hill Cimber', 'Simulated Annealing', ...
    '10 * edges', 'SA Temperature', 'K2', 'Pheno-Centric');
title('Netowrk Search Strategies: Cachexia');
ylabel('Log Likelihood Improvement, Edges, & Temperature');
xlabel('Potential States Evaluated');

%%
for i = 1:numruns
    SA_MBNet = SA_Net{i}.MakeIntoMB();
    fprintf(1,'Learning Network Parameters\n');
    SA_MBNet = LearnParamsBN(SA_MBNet);
    fprintf(1,'Predicting on Training Data\n');
    [acc, p, z] = PredictPheno(SA_MBNet, verbose);

    % output AUC or accuracy measure:
    SA_convexhullAUC(i) = AUCWorker(acc,p,z,data(:,1),doAUC);
end


toc;