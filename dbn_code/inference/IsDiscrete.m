function d = IsDiscrete(data, maxdiscvals)
%d = IsDiscrete(data, maxdiscvals)
%
% utility function that checks if the data contained in each column of DATA
% is discrete data, falling into at most MAXDISCVALS different values.
% 
% INPUT: 
% DATA: data arranged in columns
% MAXDISCVALS: maximum number of unique values a column can contain before
%   it can no longer be considered discrete.  Appropriate values are
%   dependent on the sample size, but probably between 3 and 8.
%
% OUTPUT:
% D: boolean array, true indicates that DATA(:, D(i)) is discrete.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    maxdiscvals = 10;
end

% number of columns of input:
n = size(data, 2);
d = false(1, n);

%% figure out if each variable is discrete:
r = round(data);
r = abs(r - data);
r = sum(r, 1);
for i = 1:n
    d(i) = r(i) == 0;
end
% if a node has too many values, it has to be considered continuous:
for i = 1:n
    if (d(i))
        s = length(unique(data(:,i), 'legacy'));
        if (s > maxdiscvals)
            d(i) = false;
        end
    end
end
