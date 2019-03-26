function tpval = PermTestRandomAUCWorker(data, cols, pheno)
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

phencol = strmatch(pheno, cols, 'exact');
phn = data(:,phencol);

% generate a random assignment:

z = rand(size(phn));
p = z > 0.5;
acc = sum (p == phn) / length(phn);

% turn off verbose display:
[tpval, ~] = AUCWorker(acc,p,z,phn,true,false,false);



