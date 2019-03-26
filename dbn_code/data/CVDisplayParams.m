function [m, nusigma] = CVDisplayParams(aucs, bestparams, nuVals, sigmaVals, alphaVals, maxparentsVals)
%
% take best params and try to compute average AUCs for each of the
% bestparams values
% 
% BETA
%
% copyright Michael McGeachie, 2013.

if (nargin < 3)
    % usual testing parameters for CrossValidation
    nuVals = [5, 10, 25, 50, 100];
    sigmaVals = [.3, .5, 1, 2, 3.5, 5, 10];
    alphaVals = [10, 20, 50, 100, 250];
    maxparentsVals = [2,3];
end


% find the indices for each of the best values:
% just call binincbasevec until the values match...
basevec = [length(nuVals), length(sigmaVals), length(alphaVals), length(maxparentsVals)];
inds = zeros(size(basevec));
done = false;
while (~done)
    inds = bitincbasevec(inds, basevec);
    % bitincbasevec is zero-based; just add one to make indices:
    done = ((nuVals(inds(1)+1) == bestparams.nu) && ...
           (sigmaVals(inds(2)+1) == bestparams.sigma2) && ...
           (alphaVals(inds(3)+1) == bestparams.alpha) && ...
           (maxparentsVals(inds(4)+1) == bestparams.maxParents)); 
end

% convert to actual indices:
inds = inds + 1;

% compute stepsizes array for use with ALLVALUESINDEX()
stepsizes = ones(size(basevec));
stepsizes(2:end) = cumprod(basevec(1:end-1));

m = cell(length(basevec),1);


% pick all indices where one param is locked in place:
for i = 1:length(basevec)
    m{i} = zeros(1,basevec(i));
    for j = 1:basevec(i)
        [allinds] = allvaluesindex([1:i-1,i+1:length(basevec)], basevec, stepsizes, i, j);
        %hist(aucs(allinds));
        m{i}(j) = mean(aucs(allinds));
        % also plot mean for the histogram; and mean for all AUCs, to show that
        % the mean improves...  could also compute a p-value for diff between
        % two normal distributions (although they seem to be betas).
    end
end

% also pick all indices where params (nu, sigma) are locked in place:
nusigma = zeros(basevec(1), basevec(2));
for i = 1:basevec(1)
    % max values for NU
    for j = 1:basevec(2)
        % max values for SIGMA
        [allinds] = allvaluesindex([3:length(basevec)], basevec, stepsizes, [1,2], [i,j]);
        %hist(aucs(allinds));
        nusigma(i,j) = mean(aucs(allinds));
    end
end


