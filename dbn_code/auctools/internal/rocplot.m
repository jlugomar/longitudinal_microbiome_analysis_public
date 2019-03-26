function rocplot(tp,fp)
% rocplot(tp,fp)
% plot the ROC curve with the convex hull of the ROC curve
%
% INPUT:
% tp = true positives (y values)
% fp = false positives (x values)

figure();
plot(fp,tp);

% compute the area under the ROC
r = auroc(tp,fp);
fprintf(1, 'AUROC   = %f\n', r);


idx = unique(convhull([fp;1.1], [tp;-0.1]),'legacy');
idx = idx(idx <= length(fp));
fp  = fp(idx);
tp  = tp(idx);

hold on
plot(fp, tp, 'r--');
hold off
xlabel('false positive rate');
ylabel('true positive rate');
title('ROC and ROCCH curve');

% compute the area under the ROCCH
fprintf(1, 'AUROCCH = %f\n', auroc(tp,fp));
