function bi = bitincbasevec(b,basevec)
% bi = bitincbasevec(b, basevec)
%
% takes a vector b of zero-based numbers in base BASEVEC(i), and adds one 
% to it, treating it like a multi-digit number in arbitrary mixed base.
% EXAMPLE: [1,1,0,1,0,0] -> [0,0,1,1,0,0].  (binary)
%
% returns [] when the input B represents the max value representable in
% this scheme
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.
%

if (nargin < 2)
    basevec = 2 * ones(size(b));
end

bi = b;
d = 1;
while (d <= length(b)) 
    if (bi(d) == basevec(d) -1)
        bi(d) = 0;
        d = d + 1;
    else
        break;
    end
end
if (d > length(b))
    bi = [];
else
    bi(d) = bi(d) + 1;
end