% Test DBN for Forward Simulation of most likely values:

%% load data:
version1 = '3';
fname = ['pnas_dbn_v',version1,'.txt'];
[numdata,cols] = RCSVLoad(fname,false,'\t');
data = numdata(:,2:end);
cols = cols(2:end);
subjids = numdata(:,1);


    
%% take data and split into two (overlapping) time points:
[d0,dn,nsids] = MakeTSBNData(data, subjids);

disc = IsDiscrete(data);
verbose = false;
priorPrecision.nu = 10;
priorPrecision.alpha = 10;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = 3;


%% make main call to learn the DBN:
BFTHRESH = 0;
disc = IsDiscrete(data);
searchParameter.d0 = d0;
searchParameter.dn = dn;
searchParameter.DBN = true;
searchParameter.nophenotype = true;
verbose = true;

version3 = '3';
searchParameter.annealing = false;
outfilename = ['DBN_test_v',version3];
[BNet, outstats] = FullBNLearn(data, cols, '', BFTHRESH, outfilename, ...
    priorPrecision, disc, verbose, searchParameter);
newdata = [d0,dn];
DBN3 = BNet.ConvertToDBN(newdata);


%%
fsDBN3tree = ForwardSimLearnParams(DBN3);
[predictions,likelihood,truelikelihood] = ForwardSim(DBN3, fsDBN3tree);

DBN3.WriteToGML(['DBN_testv',version3,'_DBN']);
DBN3.WriteToTGF(['DBN_testv',version3,'_DBN']);

%% get prediction and true likelihoods using PredictPartial()
preds = -1 * ones(size(DBN3.data));
z = preds;
logevprobs = preds;
trueppprobs = preds;
startcol = 5;
for i = (length(DBN3.cols)/2+1+startcol):length(DBN3.cols)
    newpheno = DBN3.cols{i};

    % make a MB:
    % need to identify a phenotype here:
    DBN3MB = DBN3.MakeIntoMB(newpheno);
    fprintf(1,'Predicting %s with a Markov Neighborhood of size %d\n', newpheno, length(DBN3MB.cols));

    % learn params
    DBN3MB = DBN3MB.LearnParams();

    % Partial Prediction:
    % enter evidence on many of the vars; predict the other vars, or output
    % their probabilities.
    match = strcmp(newpheno, DBN3MB.cols);
    predvars = find(match);
    evvars = 1:length(DBN3MB.cols);
    evvars = evvars(~match);
    [preds(:,i), z(:,i),logevprobs(:,i),trueppprobs(:,i)] = PredictPartial(DBN3MB, evvars, predvars);
    % compare = [preds,DBN4MB.GetPhenoCol()];

    % compare with ForwardSim() predictions:
    % match = strcmp(newpheno,DBN4.cols);
    % cfcomp = [preds, predictions(:,match)];
end

%% compare PredictPartial() to ForwardSim()
% or display probability of true values on bar graph.

% range of vars we predicted with both methods:
range = (length(DBN3.cols)/2+1+startcol):length(DBN3.cols);

% compare predictions for each var:
preddiffs = predictions(:,range) - preds(:,range);
numpreddiffs = sum(sum(abs(preddiffs) > 0));

% compare probabilities of inferred values:
likepreddiffs = likelihood(:,range) - z(:,range);
numlikepreddiffs = sum(sum(abs(likepreddiffs) > 0));

% compare probabilities of true data:
truepreddiffs = truelikelihood(:,range) - trueppprobs(:,range);
numtruepreddiffs = sum(sum(abs(truepreddiffs) > 0));

% result : pretty good correspondences.  Bad data for predicted likelihood
% and for "logevprobs" which seems to be uncomputable for continuous-only
% networks.

% next: plot cumulative likelihood on data graph:

save('paper_dbn_mn_predictions.mat','preds');
dlmwrite('paper_dbn_mn_predictions.txt', preds);

save('paper_dbn_predictions.mat','predictions');
dlmwrite('paper_dbn_predictions.txt', predictions);


