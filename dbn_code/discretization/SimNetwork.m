function data = SimNetwork(bn, adjmat, n)
%data = SimNetwork(bn, adjmat, n)
% 
% Take an existing bayesian network and generate data for it. Generates
% random conditional distributions for each node in the network, then
% simulates data from those distributions according to the network
% structure.
%
% INPUT:
% BN: list of NODES output by BNfromAdjMat()
% ADJMAT: adjacency matrix of the network
% N: dataset size
%
% OUTPUT:
% DATA: array of datapoints, N long by length(bn) wide.
%
%
% (c) Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.


% dataset size
if (nargin < 3)
    n = 1000;
end
data = zeros(n,length(bn));

% top sort the BN
order = TopSort(adjmat);
if (order(1) ~= 1)
    error('Top Sort did not result in first node first\n');
end

% probably need a different gaussian mean, std, for each assignment to the
% discrete nodes.

% use parameters to control the difference between those distributions
% basestd is the base standard deviation for each normal distribution
basestd = 1;
% meandiff = difference between means as a proportion of the standard
% deviation
meandiff = 1;

% count discrete nodes; starting with top sort nodes:
disc = false(1,length(order));
for i = 1:length(order)
    if (bn(order(i)).discrete)
        disc(i) = true;
    end
end

means = cell(1,length(bn));
stds = means;
params = cell(1,length(bn));

cpts = cell(1,length(order));
ndiscpvals = ones(1,length(order));
stepsize = cell(1,length(order));
discpindex = stepsize;
% go through each node
for i = 1:length(order)
    % count discrete parents of this node:
    stepsize{i} = [1];
    discpindex{i} = [];
    if (~isempty(bn(order(i)).pindex))
        for j = 1:length(bn(order(i)).pindex)
            if (bn(bn(order(i)).pindex(j)).discrete)
                discpindex{i}(end+1) = j;
            end
        end
        % for each discrete parent
        stepsize{i} = ones(1,length(discpindex{i}));
        for j = 1:length(discpindex{i})
            % check how many values that parent has
            ndiscpvals(i) = ndiscpvals(i) * length(bn(bn(order(i)).pindex(discpindex{i}(j))).values);
            if (j > 1)
                % and keep track of stepsizes:
                stepsize{i}(j) = stepsize{i}(j-1) * length(bn(bn(order(i)).pindex(discpindex{i}(j-1))).values);
            end
        end
    end
    if (bn(order(i)).discrete)
        % conditional probability table that includes probabilities for 
        % each combination of parent values.
        % this CPT can then be used to generate random data assignments
        % by stepping through the network in topological order

        % populate a CPT of that size
        cpts{i} = rand(ndiscpvals(i),length(bn(order(i)).values));
        % now normalize:
        cpts{i} = cpts{i} ./ repmat(sum(cpts{i},2),1,length(bn(order(i)).values));
    else
        % gaussian node:
        % pick a mean and std for the node, if its gaussian
        means{i} = zeros(ndiscpvals(i),1);
        stds{i} = zeros(ndiscpvals(i),1);
        ngaussp = length(bn(order(i)).pindex) - length(discpindex{i});
        params{i} = zeros(ndiscpvals(i),ngaussp);
        for j = 1:ndiscpvals(i)
            % try all nodes with same STD first:
            stds{i}(j) = basestd;
            means{i}(j) = rand(1)*6 - (j * meandiff) * stds{i}(j);
            % check number of parents, pick a parameter value for each
            % gaussian parent
            % this is an array assignment
            params{i}(j,:) = rand(1,ngaussp) * 3;
        end
    end
end

% now generate the data
% generate the head node data
r = rand(n,1);
firststep = 1 ./ length(bn(order(1)).values);
for i = 1:length(bn(order(1)).values)
    data(r >= (i-1) * firststep & r < i * firststep, bn(order(1)).colind) = bn(order(1)).values(i);
%    data(r >= .5, bn(order(1)).colind) = bn(order(1)).values(2);
%    data(r < .5, bn(order(1)).colind) = bn(order(1)).values(1);
end

% now do everything else
for i = 2:length(order)
    if (bn(order(i)).discrete)
        % discrete node
        if (isempty(bn(order(i)).pindex))
            % if there are no parents, generate random data:
            r = rand(n,1);
            base = 0;
            top = 0;
            % ok to use single-dimension indexing here with CPTS because we
            % know there are no parents, so the first dimension is 1.
            for j = 1:length(cpts{i})
                top = top + cpts{i}(j);
                data(r > base & r <= top,bn(order(i)).colind) = bn(order(i)).values(j);
                base = base + cpts{i}(j);
            end
        else
            % data has been generated for all parents
            r = rand(n,1);
            for j = 1:n
                % check values of parents
                vals = zeros(1,length(bn(order(i)).pindex));
                valinds = vals;
                for k = 1:length(bn(order(i)).pindex)
                    vals(k) = data(j,bn(bn(order(i)).pindex(k)).colind);
                    valinds(k) = find(vals(k) == bn(bn(order(i)).pindex(k)).values);
                end
                % figure out what number this is:
                ind = 1;
                for k = 1:length(valinds)
                    ind = ind + stepsize{i}(k) * (valinds(k)-1);
                end
                % finally check the CPT for this index:
                cs = cumsum(cpts{i}(ind,:));
                k = 1;
                while (r(j) > cs(k))
                    k = k + 1;
                end
                data(j,bn(order(i)).colind) = bn(order(i)).values(k);
            end
        end
    else
        % continuous node
        for d = 1:n
            % get the right discrete assignment:
            vals = zeros(1,length(discpindex{i}));
            valinds = vals;
            for k = 1:length(discpindex{i})
                vals(k) = data(d,bn(bn(order(i)).pindex(discpindex{i}(k))).colind);
                valinds(k) = find(vals(k) == bn(bn(order(i)).pindex(discpindex{i}(k))).values);
            end
            % figure out what number assignment this is:
            discassign = 1;
            for k = 1:length(valinds)
                discassign = discassign + stepsize{i}(k) * (valinds(k)-1);
            end
            if (discassign > ndiscpvals(i))
                error('somehow too many discrete values!');
            end
                
            % generate equation to generate data by looking at parents
            ind = 1;
            eq_ind = [];
            eq_coef = [];
            for j = 1: length(bn(order(i)).pindex)
                if (~bn(bn(order(i)).pindex(j)).discrete)
                    %p = params{i}(discassign,ind);
                    eq_coef(ind) = params{i}(discassign,ind);
                    eq_ind(ind) = bn(bn(order(i)).pindex(j)).colind;
                    ind = ind + 1;
                end
            end
            m = means{i}(discassign);
            st = stds{i}(discassign);
            r = random('norm',m,st);
            if (~isempty(eq_coef))
                data(d,bn(order(i)).colind) = r + eq_coef * data(d,eq_ind)';
            else
                data(d,bn(order(i)).colind) = r;
            end
        end    
    end
    
end




