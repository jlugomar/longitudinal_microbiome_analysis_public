function auc = aucmultiplot(p, z, ptrue, legendstrings, doRawROC)
%auc = aucmultiplot(p, z, ptrue, legendstrings)
%
% Makes plots of several ROCs on the same figure.  Input parameters should
% be cell arrays of ROC information, so that p{1} is the predictions for
% the first ROC, then p{2} is the second, etc...
%
% INPUT
% p : class predictions, cell of arrays.
% z : confidence of class prediction, in [0,1], cell of arrays matching p
% ptrue : true class
% legendstrings : cell array of strings to be printed in the legend of the
%   figure
% doRawROC : (optional) if true, will print both the raw ROC and the
%   convexhull of the ROC
%
%
% copyright Michael McGeachie, 2010

if (nargin < 5)
    doRawROC = false;
end

linestyles = {'--','.-','-',':'};

figure();
hold on;
for i = 1:length(p)
    % first invert the predictions on the first class:
    % (only makes sense for binary phenotype)
    zi = z{i};
    pi = p{i};
    zi(pi == 0) = 1 - zi(pi == 0);
    
    % then call the ROC routines to compute true positives, false positive rates:
    if (doRawROC)
        % raw, non-convexhull ROC
        [tp,fp] = roc(ptrue{i}, zi);
        linespec = ['b', linestyles{mod(i,length(linestyles))}];
        plot(fp, tp, linespec);
    end
    
    [tp,fp] = rocch(ptrue{i}, zi);
    auc = auroc(tp,fp);

    % plot the ROC convex hull:
    idx = unique(convhull([fp;1.1], [tp;-0.1]),'legacy');
    idx = idx(idx <= length(fp));
    fp  = fp(idx);
    tp  = tp(idx);

    linespec = ['r', linestyles{mod(i,length(linestyles))}];
    plot(fp, tp, linespec, 'LineWidth', 3);
end
% plots a line at the 50% random mark
plot([0,1],[0,1],'k:', 'LineWidth',3);
hold off;
xlabel('false positive rate (1-specificity)');
ylabel('true positive rate (sensitivity)');
legendstrings = {legendstrings{:},'random'};
legend(legendstrings,'Location','SouthEast');



