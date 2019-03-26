function mlly = LLMultiParent(xdata, ydata, alpha, ind)
% mlly = LLMultiParent(xdata, ydata, alpha, ind)
%
% Compute marginal log likelihood of the data and the structure
% xdata = data of parent
% ydata = data of child
% ind = boolean; true if parent and child are independent
%       false if child depends on parent
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.
%

% orientation
if (size(xdata,2) > size(xdata,1))
    xdata = xdata';
end
% CARDINALITY
if (~isempty(xdata))
    xvals = cell(1,size(xdata,2));
    for i = 1:size(xdata,2)
        xvals{i} = unique(xdata(:,i),'legacy');
    end
else
    xvals = [];
    xcard = 0;
end

[vinc, vind, xcard] = ValueIncrement(xvals, []);

yvals = unique(ydata,'legacy');
ycard = length(yvals);

% PRIORS
% assume priors are all uniform for now
% used to assume incoming cardinality was always: 2*3 = 6;
if (1)
    yalpha = ones(1,ycard) ./ ycard * alpha;
else
    yalpha = ones(1,ycard) .* xcard;
end

% DATA
% xdata is for the parent; ydata is for the child.
% should be same length
M = length(ydata);

% sum of alphas
sumay = sum(yalpha);

% we know this is gammaln(alpha):
galphay = gammaln(sumay);

mlly = 0;
if ind
    % y is idenpendent of x
    yMgamma = gammaln(sumay + M);
    for i = 1:ycard
        mlly = mlly + gammaln(yalpha(i) + sum(ydata == yvals(i))) - gammaln(yalpha(i));
    end
    mlly = mlly + galphay - yMgamma;
else
    % for each possible value of xdata, condition ydata on that value
    % divide the alpha parameters by each condition of x
    yalpha = yalpha ./ xcard;
    while (~isempty(vinc))
        a = repmat(vinc,size(xdata,1),1);
        matches = xdata - a;
        % matching rows are identically zero
        zs = sum(abs(matches),2);
        % select the matching rows
        ydataj = ydata(zs == 0);

        
        % divide sumalphay by # of conditions of x also:
        yMgammaj = gammaln(sumay ./ xcard + length(ydataj));
        for i = 1:ycard
            mlly = mlly + gammaln(yalpha(i) + sum(ydataj == yvals(i))) - gammaln(yalpha(i));
        end
        mlly = mlly + gammaln(sum(yalpha)) - yMgammaj;
        [vinc, vind] = ValueIncrement(xvals, vinc, vind);
    end
end
