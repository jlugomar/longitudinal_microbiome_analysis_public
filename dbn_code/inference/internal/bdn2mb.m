function mb = bdn2mb(bdnname, varname)
% MB = bdn2mb(bdnfile, mbvar) extract markov blanket 
%   from Bayesware Discoverer BDN file. 
%   varname defaults to 'phenotype'
%   '.bdn'  is appended to bdnname

if (nargin < 2)
    varname = 'phenotype';
end

net1 = bdnread(strcat(bdnname, '.bdn'));
mb = bdnmarkovblanket(net1,varname);