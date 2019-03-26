function [BN] = GetThreshBootNetwork(adjmat, data, cols, pheno, thresh, verbose)
% BN = GetThreshBootNetwork(adjmat, data, cols, pheno, thresh, verbose)
%
% function to take a continuous adjacency matrix representing the aggregate
% network generated from several bootstrap iterations, using
% BootstrapLearn.m.  Outputs an entire exhaustive network rather than just
% a Markov Blanket network.  Includes all edges that occur in greater than
% THRESH of the bootstrap networks.  Also does cycle checking and cycle
% breaking, dropping the least likely edge to break a cycle.  Outputs 
% trivial graph format files of the network if VERBOSE is TRUE.
%
% INPUT:
% ADJMAT: continuous adjacency matrix represent edge frequencies accross
%   bootstrap realizations of a CG Bayes network.
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% THRESH: the threshold for including an edge in the aggregate network.
%   Must be in [0,1], used to cut edges from non-edges in the ADJMAT. 
%   Default = 0.4, can be an array of thresholds.
% VERBOSE: will output a TGF file of the network if true
%
% OUTPUT:
% BN : a BayesNet class object representing the entire network
%
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    thresh = 0.4;
end
if (nargin < 6)
    verbose = false;
end


newmat = zeros(size(adjmat));
weightMatrix = zeros(size(adjmat));
big = 2;
while (big > thresh)
    % find highest weighted edge
    big = max(max(adjmat));
    if (big < thresh)
        break;
    end
    [mrow, mcol] = find(adjmat == big);
    for j = 1:length(mrow)        
        % and add those edges to newmat
        tempmat = newmat;
        tempwm = weightMatrix;
        tempmat(mrow(j),mcol(j)) = 1;
        tempwm(mrow(j),mcol(j)) = adjmat(mrow(j), mcol(j));
        % remove these edges from further consideration
        adjmat(mrow(j),mcol(j)) = 0;
        % then check this for cycles:
        if (hasCycle(tempmat))
            continue;
        else
            newmat = tempmat;
            weightMatrix = tempwm;
        end
    end
end

BN = BayesNet([],['ThreshNet',num2str(thresh)],newmat,weightMatrix,[],false,[],data,cols,pheno,[],{});

if (verbose)
    BN.WriteToTGF();
    outputBayesNetGraph(newmat, cols, ['Bootfile_thresh_',num2str(thresh)]);
end
