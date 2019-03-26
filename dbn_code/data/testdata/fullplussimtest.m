
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


%% full and sim anneal:
numruns = 10;

BF_THRESH = 0;
priorPrecision.maxParents = 6;
searchParameter.backtracking = true;
searchParameter.annealing = false;
searchParameter.SA_Temp_Mult = 1;
combo_lldiffs = cell(1,numruns);
combo_numedges = cell(1,numruns);
combo_numevals = cell(1,numruns);
Combo_Net = cell(1,numruns);
combo_stats1 = cell(1,numruns);
for i = 1:numruns
    fprintf(1,'Learning Network Structure w/ Exhaustive Hillclimber followed by Sim Annealing for %s\n', analysis_title);
    [Combo_Net{i}, combo_stats1{i}] = FullAndSimAnneal(data, cols, pheno, BF_THRESH, [analysis_title,'-net'], ...
        priorPrecision, disc, verbose, searchParameter);
    combo_lldiffs{i} = combo_stats1{i}.lldiffs;
    combo_numedges{i} = combo_stats1{i}.numedges;
    combo_numevals{i} = combo_stats1{i}.numevals;
end

%% plot the trajectories:
figure();
hold on;
c = colormap(jet);
[s1,~] = size(c);
ccstepsize = floor(s1 / (numruns));
h = [0 0 0];
for i = 1:numruns
    h(1) = plot(combo_lldiffs{i}, 'Color',c(ccstepsize * i,:));
    h(2) = plot(10 * combo_numedges{i},':','Color',c(ccstepsize * i,:));
end
hold off;

title('Netowrk Search Strategies: Cachexia');
ylabel('Log Likelihood Improvement, Edges, & Temperature');
xlabel('Search Space (steps taken)');


%% display as a function of number of evals
numruns = 10;
figure();
hold on;
c = colormap(jet);
[s1,~] = size(c);
ccstepsize = floor(s1 / (numruns + 4));

for i = 1:numruns
    h(1) = plot(combo_numevals{i}, combo_lldiffs{i}, 'Color',c(ccstepsize * i,:));
    h(2) = plot(combo_numevals{i}, 10 * combo_numedges{i},':','Color',c(ccstepsize * i,:));
end
hold off;


title('Netowrk Search Strategies: Cachexia');
ylabel('Log Likelihood Improvement, Edges, & Temperature');
xlabel('Potential States Evaluated');

%%
Combo_convexhullAUC = zeros(1,numruns);
for i = 1:numruns
    Combo_MBNet = Combo_Net{i}.MakeIntoMB();
    fprintf(1,'Learning Network Parameters\n');
    Combo_MBNet = LearnParamsBN(Combo_MBNet);
    fprintf(1,'Predicting on Training Data\n');
    [acc, p, z] = PredictPheno(Combo_MBNet, verbose);

    % output AUC or accuracy measure:
    Combo_convexhullAUC(i) = AUCWorker(acc,p,z,data(:,1),doAUC);
end


toc;