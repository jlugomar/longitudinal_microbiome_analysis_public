function tree = PushEvidence(tree, var, varval)
%tree = PushEvidence(tree, var, varval)
%
% Cowell's Alg 5.3: the PUSH operation (entering evidence)
% used to enter and propagate evidence of finding that variable
% VAR = value VARVAL.
%
% should initialize all tree().activeflag = true;
% and all POSTBAGS are empty.
%
% Debugging State: thoroughly tested, possible to see if exchange must be
% performed on every element in the postbag.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.



% STEP 1: enter evidence in all regressions
for i = length(tree):-1:1
    if (tree(i).index == var)
        break;
    end
    if (tree(i).activeflag)
        if (isempty(setdiff(var,tree(i).members, 'legacy')))
            for j = 1:length(tree(i).lppotential)
                tree(i).lppotential{j}(1) = tree(i).lppotential{j}(1).AddEvidence(var, varval);
            end
        end
    end
end

% STEP 2: initialize loop for pushing VARVAL to the discrete nodes
% i = elimination node for VAR, per above loop termination
n = i;
% move all potentials to the postbags:
for i = 1:length(tree(n).lppotential)
    tree(n).postbag{i} = tree(n).lppotential{i};
% maybe keeping these around won't cause problems!!
%    tree(n).lppotential{i} = [];
end
% and done with this node:
tree(n).activeflag = false;

% STEP 3: push VARVAL toward the root
while ((~isempty(tree(n).parent) && ~tree(tree(n).parent).discrete))
   for i = 1:length(tree(n).postbag)
       for j = 1:length(tree(n).postbag{i})
           % pass up to the tree(tree(n).parent) all current postbags in tree(n)
           tree(tree(n).parent) = tree(tree(n).parent).AddToPostbags(tree(n).postbag{i}(j));
       end
   end
   if (tree(tree(n).parent).activeflag)
      for i = 1:length(tree(tree(n).parent).postbag)
           % exchange:
           newlp = tree(tree(n).parent).lppotential{i}(1);
           R = tree(tree(n).parent).postbag{i}(1);
           % not sure if you have to exchange with all lps in the
           % postbag:
           [newlp, R] = Exchange(newlp, R);
           % ** Exchang Goal: R.head == VAR **
           % and secondary goal : VAR in newlp.tail
           tree(tree(n).parent).lppotential{i}(1) = newlp;
           tree(tree(n).parent).postbag{i}(1) = R;
           % enter evidence:
           tree(tree(n).parent).lppotential{i}(1) = ...
               tree(tree(n).parent).lppotential{i}(1).AddEvidence(var,varval);
       end
   end
   for i = 1:length(tree(n).postbag)
       tree(n).postbag{i} = [];
   end
   n = tree(n).parent;
end

% STEP 4: Update discrete part of tree:
% note tree(n) has a discrete parent:
for i=1:length(tree(n).postbag)
    % postbags should be length 1 here:
    % MakeWeight also includes the action of incorporating the evidence that 
    % VAR = VARVAL:
    if (isempty(tree(n).logweighttable))
        tree(n).logweighttable = zeros(1,length(tree(n).postbag));
    end
    if (i > length(tree(n).logweighttable))
        error('Overrun logweighttable in tree %d',n);
    end
    if (isempty(tree(n).postbag) || isempty(tree(n).postbag{i}))
        error('Overrun at tree %d postbag',n);
    end
    tree(n).logweighttable(i) = tree(n).logweighttable(i) + ...
        tree(n).postbag{i}(1).MakeWeight(var,varval);
    
    tree(n).postbag{i} = [];
end


