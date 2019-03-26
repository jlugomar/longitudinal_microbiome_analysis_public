function [newcpt, newcst] = factorreduceworker(cpt, varnames, evidence, cst)
%[newcpt, newcst] = factorreduceworker(cpt, varnames, evidence, cst)
% 
% factor reduce worker takes a CondProbTable and updates it with the
% evidence provided by VARNAMES being set to values EVIDENCE
%
% INPUT : 
%   CPT : CondProbTable to update
%   VARNAMES : variable names that are observed as evidence
%   EVIDENCE : values for VARNAMES, EVIDENCE(i) is a member of VARNAMES(i)
%       domain
%   CST : ClusterSetTree
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.
%

if (nargin < 4)
    cst = [];
end
newcst = cst;
newcpt = cpt;

%% for each CPT in CPTS, figure out where each varname occurs
inds = zeros(size(varnames));
for j = 1:length(varnames)
    n = varnames{j};
    for k = 1:length(cpt.factors)
        if (strcmpi(cpt.factors{k}, n))
            inds(j) = k;
            break;
        end
    end
end
% evinds is parallel index into EVIDENCE
evinds = inds > 0;
inds = inds(inds > 0);
if (isempty(inds))
    return;
end

ev = evidence(evinds);

if (~isempty(cpt.table))
    data = cpt.table(:,inds);
    a = repmat(ev,size(data,1),1);
    matches = data - a;
    % matching rows are identically zero
    zs = sum(abs(matches),2);
    % select the matching rows
    newcpt.table = cpt.table(zs == 0,:);
    % prob is now a simple number
    newcpt.logprob = cpt.logprob(zs == 0);
    if (sum(zs == 0) == 0)
        % try table-free method instead
        error('no matching table entry found for evidence');
        cpt.table  = [];
    end
    % try w/o any normalizing
    %newfs(i).prob = newfs(i).prob ./ sum(newfs(i).prob);
end
if (isempty(cpt.table))
    % compute index into values array by computational method
    % set up evidence varables index:
    evvars = cpt.factors(inds);
    % set up unassigned list of indices:
    uinds = true(size(cpt.factors));
    uinds(evvars) = false;
    unassigned = cpt.factors(uinds);
    if (isempty(cst))
        error('no ClusterSetTree given for CPT tableless index computation');
    end
    % set up stepsizes for discrete nodes
    if (isempty(cst.dvalstepsize))
        cst = cst.InitStepSizes();
    end
    % set up vallengths for discrete nodes:
    vallengths = zeros(size(cst.dvalues));
    for i = 1:length(vallengths)
        vallengths(i) = length(cst.dvalues{i});
    end
    inds = allvaluesindex(unassigned, vallengths, cst.dvalstepsize, evvars, evvarinds);
    % use these indices from the CPT.prob():
    newcpt.logprob = cpt.logprob(inds);
end

%% reduce the number of VALUES assigned to each FACTOR
for j = 1:length(inds)
    newcpt.values{inds(j)} = [ev(j)];
end

if (nargout > 1)
    newcst = cst;
end

