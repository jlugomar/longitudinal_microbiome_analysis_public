function [newpot1, newpot2] = Exchange(pot1, pot2)
%[newpot1, newpot2] = Exchange(pot1, pot2)
%
% perform Cowell's EXCHANGE operation on two lp-potentials
% EXCHANGE (pot1, pot2)
%
% INPUT : 
%   POT1 and POT2 : instantiations of class LPPotential:
%
% OUTPUT : 
%   NEWPOT1 and NEWPOT2 : results of performing exchange
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


% make sure the conditioning vars are sorted in the tails:
if (~isempty(pot1.tail))
    [pot1.tail, sorder] = sort(pot1.tail);
    pot1.params = pot1.params(sorder);
end
if (~isempty(pot2.tail))
    if (length(pot2.tail) > 1)
        [pot2.tail, sorder] = sort(pot2.tail);
        pot2.params = pot2.params(sorder);
    end       
end

% find and remove the exchange variable:
exchind = find(pot1.head == pot2.tail);
keepind = true(size(pot2.tail));
if (exchind > 0)
    keepind(exchind) = false;
    % weight previously assigned to the exchanged var:
    if (exchind > length(pot2.params))
        error('exchanging unknown index!  error in exchange!');
    end
    b = pot2.params(exchind);
    pot2.params = pot2.params(keepind);
    pot2.tail = pot2.tail(keepind);
else
    b = 0;
end

% get a list of all conditioning vars:
allvars = unique([pot1.tail, pot2.tail],'legacy');
% now extend each potential's param list with zeros for variables that do
% not condition their heads.
pot1pext = zeros(size(allvars));
pot2pext = pot1pext;
for i = 1:length(allvars)
    for j = 1:length(pot1.tail)
        if (allvars(i) == pot1.tail(j))
            pot1pext(i) = pot1.params(j);
            break;
        end
    end
    for j = 1:length(pot2.tail)
        if (allvars(i) == pot2.tail(j))
            pot2pext(i) = pot2.params(j);
            break;
        end
    end
end

newpot1 = pot1;
newpot2 = pot2;

% set up the simple (second) potential:
newpot2.params = pot2pext + b * pot1pext;
newpot2.const = pot2.const + b * pot1.const;
newpot2.sigma = pot2.sigma + b^2 * pot1.sigma;

% other potential is trickier, based on the values of SIGMA.
% first add the new parameterization on the exchanged var:
if (pot2.sigma == 0 && pot1.sigma == 0)
    % in this case newpot1.params doesn't change, for that matter, neither
    % does sigma:
    newpot1.sigma = 0;
    % in this case, doesn't depend on the new exchanged var:
    newpot1.params(end + 1) = 0;
elseif (newpot2.sigma == 0)
    newpot1.sigma = 0;
    newpot1.params = -1 * pot2pext / b;
    newpot1.params(end + 1) = 1/b;
    newpot1.const = -1 * pot2.const / b;
else
    % this is the most general/common case
    k = pot2.sigma + b^2 * pot1.sigma;
    newpot1.sigma = pot1.sigma * pot2.sigma / k;
    newpot1.params = (pot1pext * pot2.sigma - b * pot2pext * ...
        pot1.sigma) / k;
    newpot1.params(end + 1) = b * pot1.sigma / k;
    newpot1.const = (pot1.const * pot2.sigma - b * pot2.const * ...
        pot1.sigma) / k;
end

newpot1.tail = [allvars, pot2.head];
newpot1.tail = newpot1.tail(newpot1.params ~= 0);
newpot1.params = newpot1.params(newpot1.params ~= 0);

newpot2.tail = allvars;
newpot2.tail = newpot2.tail(newpot2.params ~= 0);
newpot2.params = newpot2.params(newpot2.params ~= 0);





