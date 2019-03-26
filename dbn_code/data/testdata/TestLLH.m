% test loglikelihood:

% test on chachexia dataset:
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

% first log xform all metabolites in data2:
FIRSTMET = 2; % looks like first metabolite is in the 2nd column
for i = FIRSTMET:size(data,2)
    data(:,i) = log2(data(:,i));
end


%% check loglikelihood, compare vs. unnormalized
analysis_title = 'Likelihood Comparison of Chachexia';
fprintf(1,'\n\n======  start new analysis: %s  =====\n\n',analysis_title);
tic;

% common parameter values:
priorPrecision.nu = 10;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 6;

LLH1 = EmptyLogLikelihood(data, disc, priorPrecision);

% now normalize:
disc = IsDiscrete(data);
for i = 1:length(disc)
    if (~disc(i))
        % normalize:
        md = mean(data(:,i));
        stdd = std(data(:,i));
        data(:,i) = (data(:,i) - md) ./ stdd;
    end
end

LLH2 = EmptyLogLikelihood(data, disc, priorPrecision);

toc;
fprintf(1,'\nUnnormalized LLH: %f, Normalized LLH: %f\n\n',LLH1,LLH2);


