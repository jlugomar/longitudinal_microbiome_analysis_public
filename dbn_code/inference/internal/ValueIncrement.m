function [vinc, vind, numvals] = ValueIncrement(vals, v, vind)
%[vinc, vind, numvals] = ValueIncrement(vals, v, vind)
%
% takes a cell of arrays of numeric (discrete) values
% and an array that is an assignment to each member of VALS
% so that V(i) is a member of VALS{i}
% 
% returns the value array that is obtained from V() by incrementing V(i) if
% i is the first element of V(i) that can be incremented w/o going out of
% bounds on VALS{i}.
%
% this does an "ALL COMBINATIONS" iteration through assignments to VALS{}
% with successive calls to VALUEINCREMENT() like:
% v = ValueIncrement(vals, v);
%
% if input V = [], will return the first assignment to VALS{}, and will
% assign this output var:
% output : NUMVALS = total number of different value assignments possible
%
% Debugging Status: TESTED SEPARATELY.

if (nargin < 3)
    vind = ones(size(v));
end
numvals = 1;
if (isempty(v))
    vinc = zeros(size(vals));
    for i = 1:length(vals)
        numvals = length(vals{i}) * numvals;
        v1 = vals{i};
        vind(i) = 1;
        if (iscell(v1))
            vinc(i) = v1{1};
        else
            vinc(i) = v1(1);
        end
    end
    return;
end
vinc = v;
if (nargin < 3)
    % set up the indices
    for i = 1:length(v)
        for j = 1:length(vals{i})
            % how to do this as cells?
            if (vals{i}{j} == v(i))
                vind(i) = j;
                break
            end
        end
    end
end
% now do the increment:
i = 1;
while i <= length(v)
    if (vind(i) == length(vals{i}))
        vind(i) = 1;
        v1 = vals{i};
        if (iscell(v1))
            vinc(i) = v1{1};
        else
            vinc(i) = v1(1);
        end
        i = i + 1;
    else
        break;
    end
end
if (i > length(v))
    vinc = [];
else
    v1 = vals{i};
    vind(i) = vind(i) + 1;
    if (iscell(v1))
        vinc(i) = v1{vind(i)};
    else
        vinc(i) = v1(vind(i));
    end
end
