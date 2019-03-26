function p = SNP_HWE(obs_hets, obs_hom1, obs_hom2)
% p = SNP_HWE(obs_hets, obs_hom1, obs_hom2)
%
% This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
% Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
% Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  
%
% NOTE: return code of -1.0 signals an error condition
%
% translation to MATLAB by mike mcgeachie july 2008.

if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
    p = -1;
    return;
end

% total number of genotypes
N = obs_hom1 + obs_hom2 + obs_hets;
if (N <= 0)
    p = -1;
    return;
end

% rare homozygotes, common homozygotes
obs_homr = min(obs_hom1, obs_hom2);
obs_homc = max(obs_hom1, obs_hom2);

% number of rare allele copies
rare = obs_homr * 2 + obs_hets;

% Initialize probability array
probs = zeros(1, 1 + rare);

% Find midpoint of the distribution
mid = floor(rare * ( 2 * N - rare) / (2 * N));
if ( mod(mid,2) ~= mod(rare, 2) ) 
    mid = mid + 1;
end

probs(mid + 1) = 1;
mysum = 1;

% Calculate probablities from midpoint down 
curr_hets = mid;
curr_homr = (rare - mid) / 2;
curr_homc = N - curr_hets - curr_homr;

while ( curr_hets >=  2)
  probs(curr_hets - 1) = probs(curr_hets + 1) * curr_hets * (curr_hets - 1.0) ...
      / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
  mysum = mysum + probs(curr_hets - 1);

  % 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
  curr_hets = curr_hets - 2;
  curr_homr = curr_homr + 1;
  curr_homc = curr_homc + 1;
end

% Calculate probabilities from midpoint up
curr_hets = mid;
curr_homr = (rare - mid) / 2;
curr_homc = N - curr_hets - curr_homr;

while ( curr_hets <= rare - 2)
  probs(curr_hets + 3) = probs(curr_hets + 1) * 4.0 * curr_homr * curr_homc ...
      / ((curr_hets + 2.0) * (curr_hets + 1.0));
  mysum = mysum + probs(curr_hets + 3);

  % add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
  curr_hets = curr_hets + 2;
  curr_homr = curr_homr - 1;
  curr_homc = curr_homc - 1;
end

% P-value calculation
target = probs(obs_hets + 1);
p = min(1, sum(probs(probs <= target))/ mysum);

