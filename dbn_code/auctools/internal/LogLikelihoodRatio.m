function LLH = LogLikelihoodRatio(preds,class)
%LLH = LogLikelihoodRatio(preds,class)
%
% Computes loglikelihood ratio for predictions on a binary class.
% Take predictions and true class in and compute the likelihood of this
% data according to the prediciton probabilities
% preds must be % confidences between 0 and 1
%  

% transform preds by inverting zero predictors
preds(class == 0) = 1 - preds(class == 0);

% log likelihood of the data given the model
LLHm = sum(log(preds));

% log likelihood of teh data given a null-model:
freq = sum(class)/length(class);
% just say that all classes are equally likely, apriori
freq = .5;

% the null model just guesses at the class with probability = freq

prednull = ones(size(preds)) * freq;
prednull(class == 0) = 1 - prednull(class == 0);
LLHnull = sum(log(prednull));

LLH = LLHm / LLHnull;
