function ll = loglikelogreg(x,y)

% first fit the model:
warning off;
[b,~,~] = glmfit(x,y,'binomial','link','logit');
warning on;

% predict on the data
yzs = glmval(b,x,'logit');
if (size(yzs,2) > 1)
    yzs = yzs(:,2);
end

% make predictions : 
ypred = yzs > .5;

% flip probs on class = zero guesses and on incorrect cases:
yprobs = yzs;
yprobs(ypred == 0) = 1-yprobs(ypred == 0);
yprobs(ypred ~= y) = 1-yprobs(ypred ~= y);

ll = log(yprobs);
ll = sum(ll);