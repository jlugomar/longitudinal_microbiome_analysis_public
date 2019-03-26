function outdata = SimulateData(numPatients, numSNPs, pheno, skipfileout, skipmlltest)
%
% make fake data with the Minor Allele Frequencies similar to our observed
% MAFs
%
% for now, use random draws from Beta(3,9) to simulate MAF, and take that
% result MOD 0.5.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    skipmlltest = 0;
end    
if (nargin < 4)
    skipfileout = 0;
end
if (nargin < 3 || isempty(pheno))
    skippheno = true;
    skipmlltest = true;
    skipfileout = true;
end


mafs = random('beta',3,9,[1,numSNPs]);
mafs(mafs > 0.5) = mafs(mafs > 0.5) - 0.5;
fakedata = zeros(numPatients, numSNPs);

% now generate numPatients draws from that MAF for each chromosome:
for i = 1:numSNPs
    r = rand(numPatients,1);
    q2 = mafs(i)^2;
    p2 = (1-mafs(i))^2;
    hom1 = (r >= (1-p2));
    hom2 = (r <= q2);
    hets = (r > q2 & r < (1-p2));
    fakedata(hom1,i) = 1;
    fakedata(hets,i) = 2;
    fakedata(hom2,i) = 3;
end

% do filters for MAF < 5% and HWE failure:

goodmaf = ~MAFCheck(fakedata);

sim_hets = sum(fakedata == 2);
sim_hom1 = sum(fakedata == 1);
sim_hom2 = sum(fakedata == 3);

hwe = zeros(1,numSNPs);
for i = 1:numSNPs
    hwe(i) = SNP_HWE(sim_hets(i), sim_hom1(i), sim_hom2(i));
end

fakedata = fakedata(:, (goodmaf & (hwe > 0.05)));
if (~skippheno)
    fakedata = [pheno, fakedata];
end
if (skipmlltest)
    outdata = fakedata;
else    
    mll = CheckMLL(fakedata);
    inds = find (mll > 0);
    [sortmll, maxedinds] = sort(mll(inds),2,'descend');
    inds = inds(maxedinds);
    outdata = [pheno, fakedata(:,inds)];
end
[s1,s2] = size(outdata);

if (skipfileout ~= 0)
    return;
end

titlearray = cell(s2,1);
titlearray{1} = 'phenotype';
for i=2:s2
    titlearray{i} = ['bs' num2str(i)];
end

fileoutworker('simdata3.txt', titlearray, outdata);
