function aucs = SimLogFit(data, cols, pheno, snpnames, numsims)
%aucs = SimLogFit(data, cols, pheno, snpnames, numsims)
%
% perform linear regression model building on many iterations of simulated
% SNP data.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

% find pheno col
phencol = strmatch(pheno, cols, 'exact');
xind = true(size(cols));
xind(phencol) = false;
y = data(:,phencol);

% find SNP columns:
snpcols = zeros(size(snpnames));
for i = 1:length(snpnames)
    snpcols(i) = strmatch(snpnames{i}, cols, 'exact');
end

aucs = zeros(1,numsims);
for i = 1:numsims
    successsimulated = 0;
    while (successsimulated ~= sum(snpcols > 0))
        fakesnps = SimulateData(length(y), sum(snpcols > 0));
        successsimulated = size(fakesnps,2);
    end
    fakesnps = fakesnps - 1;

    fdata = data;
    fdata(:, snpcols(snpcols > 0)) = fakesnps;
    fx = fdata(:,xind);
    warning off;
    [b,~,stats] = glmfit(fx,y,'binomial','link','logit');
    warning on;
    % use the values that are statistically significant:
    good = stats.p(2:end) <= 0.05;
    if (stats.p(1) <= 0.05)
        const = 'on';
    else
        const = 'off';
    end

    if (sum(good) > 0)
        % predict on the data
        yzs = glmval(b(stats.p <= 0.05),fx(:,good),'logit','constant',const);
        if (size(yzs,2) > 1)
            yzs = yzs(:,2);
        end
    else
        yzs = zeros(size(y));
    end
    ypred = yzs > .5;
    acc = ypred == data(:,phencol);
    % compute AUC of these predictions:
    aucs(i) = AUCWorker(acc,ypred,yzs,fdata(:,phencol),true,true,false);
end

