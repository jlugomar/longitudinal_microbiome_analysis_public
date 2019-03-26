function tpval = adaptivepermtester(numsims, testval, testdir, fhandle, varargin)
%tpval = adaptivepermtester(numsims, testval, testdir, fhandle, varargin)
%
% This function computes p-values for a teststatistic using permutation
% testing for a maximum of numsims trials.  It will abandon perm testings
% for values that do not look promising.
%
% INPUT:
% NUMSIMS: maximum number of simulations to perform for p-value computation
% TESTVAL: the value of the test on the real data
% TESTDIR: if >0, means to do a "greater than" test on TESTVAL > X;
%   otherwise to do a "less than" test on TESTVAL < X
% FHANDLE: a funciton handle containing the name of a function to call with input
% arguments VARARGIN.  use @fname to generate function handles.
%
% OUTPUT: 
% TPVAL: p-value of the TESTVAL relative to the simulated null distribution
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

killthresh = .5;
iterations = 10;
simvals = zeros(1,numsims);
i = 0;
while i < numsims
    i = i + 1;
    % do one permutation:
    simvals(i) = fhandle(varargin{:});
    if (i == iterations)
        if (testdir)
            s = sum(simvals(1:i) > testval) / i;
        else
            s = sum(simvals(1:i) < testval) / i;
        end
        if (s > killthresh)
            break;
        else
            iterations = iterations * 2;
            if (iterations > numsims)
                iterations = numsims;
            end
            killthresh = .1 + .4 ^ (1 + 0.01 * iterations);
        end
    end
end
if (testdir)
    tpval = sum(simvals(1:i) > testval) / i;
else
    tpval = sum(simvals(1:i) < testval) / i;
end
