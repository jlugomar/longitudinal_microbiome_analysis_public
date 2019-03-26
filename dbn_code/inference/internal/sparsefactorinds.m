function inds = sparsefactorinds(f, varind, varval)
%
% rather than holding a whole table of combinatorial index comutations :
%   (ie, it looks like : 
%       0 0 0 : x1
%       0 0 1 : x2
%       0 1 0 : x3, ... etc )
% this computes the appripriate indices for variables that are assigned
% values varval.
%
% INPUT : 
%   F : CondProbTable
%   VARIND : indices into F.factors{VARIND}
%   VARVAL : value of f.factors{varind} to match VARVAL
%
% returns an index vector into f.probs for values of f.factors{varind}
% matching VARVAL
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


vlengths = zeros(size(f.values));
for i = 1:length(f.values)
    vlengths(i) = length(f.values{i});
end
cp = cumprod(vlengths);

vrank = 0;
for i = 1:length(f.values{varind})
    if (f.values{varind}(i) == varval)
        vrank = i;
        break;
    end
end

if (varind > 1)
    fbstart = (vrank-1) * cp(varind-1) + 1;
    fbend = vrank * cp(varind-1);
    blockend = cp(varind);
    inds = false(1,blockend);
    inds(fbstart:fbend) = true;
    
    % slickness: 
    repeats = cp(end) / cp(varind);
    inds = repmat(inds,1,repeats);
else
    inds = false(1,vlengths(varind));
    inds(vrank) = true;
    repeats = cp(end) / cp(varind);
    inds = repmat(inds,1,repeats);
end
    
    
    
    
    
    
    
    
    
    
    
    
    