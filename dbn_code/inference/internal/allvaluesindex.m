function [inds] = allvaluesindex(unassigned, vallengths, valstepsizes, evvars, evvalinds)
%[inds] = allvaluesindex(unassigned, vallengths, valstepsizes, evvars, evvalinds)
%
% One in a series of files that deals with assigning indices to numbers in
% arbitrary mixed base representations.  This is useful for indexing
% through an orderd set of arrays, each of a different length.
% This function returns ALL INDICES that have the EVVARS fixed at values
% indexed by EVVALINDS.
%
% UNASSIGNED is an array of indices into the VALLENGTHS and VALSTEPS 
%   corresponding to discrete vars that are not observed.  The output will
%   be indices that indicate all values for the UNASSIGNED.
% VALLENGTHS(UNASSIGNED(i)) is the # of possible values UNASSIGNED(i) may obtain
% VALSTEPSIZES(UNASSIGNED(i)) is the stepsize associated with UNASSIGNED(i) 
%   ~ similar to the cumulative product of VALLENGTHS.
%   actually = [1, STEPS(1:end)] where:
%       STEPS = cumprod(VALLENGTHS(1:end-1))
% EVVARS(i) means that MEMBER(EVVAR(i)) has been observed to have a 
%   particular fixed value.
% EVVALINDS(i) is the index of the value observed for MEMBER(EVVAR(i))
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% error checking:
if (isempty(evvars))
    error('No Observations for Evvars!');
elseif (round(evvars(1)) ~= evvars(1))
    error('Discrete values expected!');
elseif (evvars(1) == 0)
    warning('Zero Index!');
end
% more error checking:    
ers = evvalinds > vallengths(evvars);
if (sum(ers) > 0)
    valmaxes = vallengths(evvars);
    evvalinds(ers) = valmaxes(ers);
end
increment = valstepsizes(evvars) * (evvalinds - 1)';
cp = cumprod(vallengths(unassigned));
if (isempty(cp))
    inds = [];
else
    inds = zeros(cp(end),1);
end
% loop through assignments to (unassigned)
i = 1;
if (~isempty(unassigned))
    bvec = zeros(size(unassigned));
    while (~isempty(bvec))
        inds(i) = bvec * (valstepsizes(unassigned))';
        % ind is computed for a zero-based array, so requies a +1 adjustment
        inds(i) = inds(i) + increment + 1; 
        bvec = bitincbasevec(bvec,vallengths(unassigned));
        i = i + 1;
    end
else
    % when unassigned is empty, all that remains is the increment : 
    inds = increment + 1;
end
