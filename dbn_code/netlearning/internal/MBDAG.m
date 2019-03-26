function newmat = MBDAG(contadjmat, phncol)
%newmat = MBDAG(contadjmat, phncol)
%
% Take a continuous adjacency matrix and use it to make a DAG by breaking
% cycles by dropping the least weighted edge in the CONTADJMAT
%

% start with phncol:
enum = 1:length(contadjmat);
children = enum(contadjmat(phncol,:) > 0);

newmat = zeros(length(contadjmat));
newmat(phncol,:) = contadjmat(phncol,:);
% pick a child and check its parents:
for i = 1:length(children)
    c = children(i);
    parents = enum(contadjmat(:,c) > 0);
    for j = 1:length(parents)
        p = parents(j);
        newmat(p,c) = contadjmat(p,c);
        if (hasCycle(newmat))
            newmat(p,c) = 0;
        end

        %check each child of this parent to see if its 
        pcs = enum(contadjmat(p,:) > 0);
        %skip c
        for k = 1:length(pcs)
            pchild = pcs(k);
            newmat(p, pchild) = contadjmat(p,pchild);
            if (hasCycle(newmat))
                newmat(p, pchild) = 0;
            end
        end           
    end
end

