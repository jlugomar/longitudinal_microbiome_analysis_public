function b = NaiveBayes2(data, ind)
% b = NaiveBayes2(data, ind)
%
% return a list of assignments for the given data(:,ind) tested (and
% trained) on data(:,1).  Assignments are probability values between [0,1].
%
% Uses a cheap, quick approximation to niave bayes; no bayesian priors and
% only works with discrete variables.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

% constant:
MAXVALS = 20;

% output vector
b = zeros(size(data(:,1)));

% figure out number of values of the attribute
valinds = zeros(MAXVALS,1);
for j = 0:MAXVALS-1
    a = sum(data(:,ind) == j);
    valinds(j+1) = a > 0;
end
vals = 0:MAXVALS-1;
vals = vals(logical(valinds));

% compute conditional distribution for each value

baserate = sum(data(:,1))/length(data(:,1));
guess = zeros(size(vals));
for j = 1:length(vals);
    totalexamples = sum(data(:,ind) == vals(j));
    posexamples = sum(data(data(:,ind) == vals(j),1));
    posprob = baserate * posexamples / totalexamples;
    negexamples = sum(data(data(:,ind) == vals(j),1) == 0);
    negprob = baserate * negexamples / totalexamples;
    
    %scale up probabilities so they sum to one:
    probsum = posprob + negprob;
    posprob = posprob / probsum;
    negprob = negprob / probsum;
    guess(j) = max(posprob,negprob);
end

% now just loop (stupidly) through every data point to get accuracy
% of the classifier
for k = 1:length(data(:,ind))
    index = find(vals == data(k,ind)); % error here is with value == -1
    if (isempty(index))
        % guess here?
        % no information, so use the base rate:
        b(k) = baserate;
        continue;
    end
    b(k) = guess(index);
end


