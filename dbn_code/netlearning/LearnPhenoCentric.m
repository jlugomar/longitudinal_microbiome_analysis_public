function [MBNet, FullNet, outstats]=LearnPhenoCentric(data, cols, pheno, priorPrecision, BF_THRESH, verbose, disc)
% [MBNet, FullNet, outstats]=LearnPhenoCentric(data, cols, pheno, priorPrecision, BF_THRESH, verbose, disc)
%
% learns a Pheno-centric BN that is based around the phenotype.  All nodes
% are checked to see if the phenotype is a good parent, and then those
% nodes are checked for additional parents.
%
% INPUT:
% DATA: data array (data points X by variables)
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array.
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% BF_THRESH: log Bayes Factor threshold for edge inclusion.  Edges below this
%   threshold will not be created.  This is the main stopping criteria.
%   Default = 0.
% VERBOSE: if true, increases output.
% DISC: (optional) can be used to specify the columns that contain discrete
%   data.
%
% OUTPUT:
% MBNet: Class BayesNet object that represents the markov blanket bayesian
%   network
% FullNet: Class BayesNet object that represents the full bayesian network
% OUTSTATS: structure with fields describing the characteristics of the
%       search procedure, in arrays per "step;" a step is either an edge
%       added or removed.
%   outstats.lldiffs: difference in loglikelihood at each step of the algorithm.
%   outstats.numedges: number of edges in the network at each step of the algorithm.
%   outstats.numevals: number of potential network states evaluated at each
%       step.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 6)
    verbose = false;
end
if (nargin < 7)
    disc = IsDiscrete(data);
end
if (nargin < 5)
    BF_THRESH = 0;
end

% find pheno columhn:
phncol = find(strcmp(pheno, cols));
phendata = data(:,phncol);
adjmat = zeros(length(cols));
weightMatrix = zeros(length(cols));

%% learn nodes for which pheno is a good parent
if (verbose)
    disp('Learn nodes with pheno as parent!');
end
tempLogLLH_Ind = zeros(1,length(cols));
tempLogLLH_Dep = zeros(1,length(cols));
numevals = 0;
nedges = 0;
evals = 0;
lldiffs = 0;
addededge = {};
for v=1:length(cols)
    if (v == phncol)
        continue;
    end
    if (disc(v))
        if (disc(phncol))
            [tempModel]=learnLocalBN_DiscToDisc([],data(:,v)',priorPrecision);
            tempLogLLH_Ind(v) = tempModel.logLLH;
            [tempModel]=learnLocalBN_DiscToDisc(phendata',data(:,v)',priorPrecision);
            tempLogLLH_Dep(v) = tempModel.logLLH;    
        else
            % can't have a disc child with continuous parent, even if that
            % parent is the phenotype
            continue;
        end
    else
        if (disc(phncol))
            [tempModel]=learnLocalBN_DiscToCont([], data(:,v)', priorPrecision);
            tempLogLLH_Ind(v) = tempModel.logLLH;
            [tempModel]=learnLocalBN_DiscToCont(phendata', data(:,v)', priorPrecision);    
            tempLogLLH_Dep(v) = tempModel.logLLH;    
        else
            [tempModel]=learnLocalBN_ContToCont([], data(:,v)', priorPrecision);
            tempLogLLH_Ind(v) = tempModel.logLLH;
            [tempModel]=learnLocalBN_ContToCont(phendata', data(:,v)', priorPrecision);    
            tempLogLLH_Dep(v) = tempModel.logLLH;    
        end            
    end
    numevals = numevals + 1;
    
    % if Dep > Ind, add edge (pheno, v) to the network
    if (tempLogLLH_Dep(v)-tempLogLLH_Ind(v) > BF_THRESH)
        adjmat(phncol, v) = 1;
        weightMatrix(phncol, v) = tempLogLLH_Dep(v)-tempLogLLH_Ind(v);
        nedges(end+1) = nedges(end) + 1;
        evals(end+1) = numevals;
        lldiffs(end+1) = lldiffs(end) + tempLogLLH_Dep(v)-tempLogLLH_Ind(v);
        addededge{end+1} = [phncol,v];
    end
end

%% now look through all these nodes and see if they need other parents
for v = 1:length(cols)
    if (v == phncol)
        continue;
    end
    if (adjmat(phncol, v))
        done = false;
        while(~done && sum(adjmat(:,v)) < priorPrecision.maxParents)
            llh = -Inf * ones(1,length(cols));
            for p = 1:length(cols)
                % get info on current parents
                curp = adjmat(:,v)';
                curp = logical(curp);
                curpdisc = disc(curp);
                curpdata = data(:,curp);

                % if this parent is already a parent, don't consider it
                if (curp(p))
                    continue;
                end
                % can't be your own parent:
                if (p == v)
                    continue;
                end
                % see if p is a good parent for v
                if (disc(v) && ~disc(p))
                    continue;
                end

                % do a comparison 
                if (disc(v))
                    [baseModel]=learnLocalBN_DiscToDisc(curpdata', data(:,v)', priorPrecision);
                    [pModel]=learnLocalBN_DiscToDisc([curpdata, data(:,p)]', data(:,v)', priorPrecision);
                else
                    [baseModel]=learnLocalBN_MixToCont(curpdata(:,~curpdisc)',...
                        curpdata(:,curpdisc)', data(:,v)', priorPrecision);
                    if (disc(p))
                        [pModel]=learnLocalBN_MixToCont(curpdata(:,~curpdisc)', ...
                            [data(:,p), curpdata(:,curpdisc)]', data(:,v)', priorPrecision);
                    else
                        [pModel]=learnLocalBN_MixToCont(...
                            [curpdata(:,~curpdisc), data(:,p)]',curpdata(:,curpdisc)',data(:,v)',priorPrecision);
                    end
                end
                numevals = numevals + 1;
                llh(p) = pModel.logLLH - baseModel.logLLH;
            end %for
            
            % check max value
            [best, bind] = max(llh);
            if (best > BF_THRESH)
                tempadjmat = adjmat;
                tempWeightMat = weightMatrix;
                tempadjmat(bind, v) = 1;
                tempWeightMat(bind,v) = best;
                % unless we created a cycle
                iscycle = dfsCycleCheck(tempadjmat,v);
                while (iscycle)
                    llh(bind) = -Inf;
                    [best,bind] = max(llh);
                    if (best <= BF_THRESH)
                        done = true;
                        iscycle = true;
                        break;
                    end
                    tempadjmat = adjmat;
                    tempWeightMat = weightMatrix;
                    tempadjmat(bind, v) = 1;                    
                    tempWeightMat(bind,v) = best;
                    iscycle = dfsCycleCheck(tempadjmat,v);
                end
                if (~iscycle)
                    % really add this parent to the adjmat:
                    adjmat = tempadjmat;
                    weightMatrix = tempWeightMat;
                    nedges(end+1) = nedges(end) + 1;
                    lldiffs(end+1) = best + lldiffs(end);
                    evals(end+1) = numevals;
                    addededge{end+1} = [bind,v];
                end
            else
                % end the search
                done = true;
            end
        end %while
    end
end
% end learning network 


%% output, using new BayesNet class object:
FullNet = BayesNet([],'',adjmat,weightMatrix,[],true,disc,data,cols,pheno, priorPrecision,{});
% have to invoke MBNet.MakeIntoMB() to limit the domain:
MBNet = FullNet.MakeIntoMB(); 

outstats.numevals = evals;
outstats.lldiffs = lldiffs;
outstats.numedges = nedges;
outstats.addededge = addededge;

return;