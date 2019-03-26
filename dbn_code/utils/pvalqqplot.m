function [gif, cout, sorder] = pvalqqplot(pvals, logs, colorscale, nogif)
% [gif] = pvalqqplot(pvals, logs, docorrect)
% 
% do a log-log QQ plot of p-values vs. a uniform distribution
%
% INPUT:
% PVALS: list of pvalues
% LOGS: boolean, should be "true", plots negative log pvalues.  Can leave
%   blank.
% COLORSCALE: will plot x's based on color specified by this value;
%   parallel to PVALS
% NOGIF: boolean, by defaut FALSE.  If true, will not print the GIF and
%   LAMBDA on the graph.
%
% OUTPUT: 
% GIF: Genome-Inflation Factor.  This is a measure of the stratification
%   within the population; or the general systematic enrichment of the
%   p-values for significance.  If (much) greater than 1.0, it means the null
%   distribution does not hold; and therefor the modeling assumptions are
%   incorrect.
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    logs = true;
end
docolor = true;
if (nargin < 3)
    docolor = false;
end
if (nargin < 4)
    nogif = false;
end

% flip input
if (size(pvals,1) ~= 1)
    pvals = pvals';
end

if (logs)
    sample = (1:length(pvals)) ./ (length(pvals) +1);
    x = -1*log10(sample);
    [y,sorder] = sort(pvals);
    y = -1*log10(y);
else
    step = (max(pvals) - min(pvals))/1000;
    sample = min(pvals):step:max(pvals);
    x = sample;
    [y,sorder] = sort(pvals);
end
clear sample;
figure();
hold on;
grid;
%maxh=ceil(log10(length(x)));
if (docolor)
    colormap('jet');
    cout = colorscale(sorder);
    scatter(x,y,50,cout);
else
    plot(x,y,'bx');
end
plot(x,x,'r-');

% compute genomic inflation factor:
gif = (x * y') / (x * x');
xlabel('Expected -log10 p-values');
ylabel('Observed -log10 p-values');
if (~nogif)
    text(max(x)/50, max(y), ['  GIF = ',num2str(gif)], 'FontSize',14);
end
% compute lambda:
lambda = lambdaGIF(pvals);
if (~nogif)
    text(max(x)/50, max(y) - max(y)/12, ['  lambda = ',num2str(lambda)], 'FontSize',14);
end
clear pvals;


