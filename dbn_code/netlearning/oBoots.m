function adjmat = oBoots(fname, nboots, priorPrecision, ALG, disc)
% Do Bootstrapping for Metabolomics on Orchestra, or other parallel cluster.

%run ('bnpathscript.m');

if (nargin < 2)
    nboots = 1;
end
if (nargin < 3)
    priorPrecision.nu = 10;
    priorPrecision.sigma2 = 1;
    priorPrecision.alpha = 10;
    priorPrecision.maxParents = 4;
end
if (nargin < 4)
    ALG = 3; % exhaustive search
end
    


% seed random number generator here so that we get different sequencces of
% numbers from each instantiation of this on parallel machines:
% Aparently the issue that each matlab worker process initializes the
% random number generator by seeding on it's job id, which is just an
% integer from 1-n.  So multiple runs of the same parallel job result in
% the same sequences of random numbers.
rng('shuffle');

% load data:
[data, cols] = RCSVLoad(fname,true);
pheno = cols{1};
if (nargin < 5)
    disc = IsDiscrete(data);
end
BFTHRESH = 0;

task = getCurrentTask();
errorfile = ['task_',task.Name,'_error.txt'];
fid = fopen(errorfile,'w');
fprintf(fid, 'Reading main datafile: %s ... \n', fname);
fprintf(fid, 'Loaded data: size %d by %d\n', size(data,1), size(data,2));
fprintf(fid, 'Loaded cols: size %d, with pheno : %s\n', length(cols), pheno);
fprintf(fid, 'doing N = %d bootstraps.\n', nboots);
fprintf(fid, '\tpriorPrecision.nu = %d\n', priorPrecision.nu);
fprintf(fid, '\tpriorPrecision.sigma2 = %d\n', priorPrecision.sigma2);
fprintf(fid, '\tpriorPrecision.alpha = %d\n', priorPrecision.alpha);
fprintf(fid, '\tpriorPrecision.maxParents = %d\n', priorPrecision.maxParents);
fclose(fid);

% and do the bootstrapping:
adjmat = BootstrapLearn(data, cols, pheno, priorPrecision, nboots, ALG, false,...
    BFTHRESH, disc);

