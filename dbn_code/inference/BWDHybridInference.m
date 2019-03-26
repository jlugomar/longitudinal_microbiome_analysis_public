function [acc,p] = BWDHybridInference(bwfilename, datafile, predvar, mbfile, priorPrecision)
%[acc,p] = BWDHybridInference(bwfilename, datafile, predvar, mbfile, priorPrecision)
%
% read a BayesWare Discoverer network file into MATLAB and do inference on
% the datafile using the hybrid network algorithm.
%
% INPUT:
% BWFILENAME: filename for the BayeswareDiscoverer network to read in
% DATAFILE: filename of the data matching BWFILENAME to predict upon
% PREDVAR: name of phenotype variable to predict
% MBFILE: can be false; will generate a file of just the markov blanket of
%   the PREDVAR of the BWFILENAME network.
% PRIORPRECISION: structure where:
%   priorPrecision.alphs: should be set; prior frequency of discrete
%   evidence.
%
% OUTPUT:
% ACC: accuracy of predictor
% P: predictions per instance of the PREDVAR variable
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    priorPrecision.alpha = 6;
end
if (nargin < 4)
    mbfile = false;
end
if (nargin < 3)
    predvar = 'PHENOTYPE';
end

% if we weren't given a Markov-Blanket file, create one first:
if (~mbfile)
    fprintf(1,'Reading data from BayesWare Discoverer file %s\n', bwfilename);
    % simplify (the same) Bayes Net to a markov blanket
    net1 = bdnread(bwfilename);
    mb = bdnmarkovblanket(net1,predvar);

    % save a copy of this MB-bayesnet in a BDN file
    bwfilename = bdnmbwrite(bwfilename, mb, predvar);
end

% read the markov-blanketed BWD file:
mbnodes = bdnread(bwfilename);

fprintf(1,'Reading data from %s\n', datafile);
[data,cols] = ReadHeaderDataFile(datafile);
fprintf(1,'Learning Network Parameters\n');
BN = BayesNet(mbnodes, '', [],[],[],[],[],data,cols,predvar,priorPrecision);
BN = LearnParams(BN);
fprintf(1,'Predicting on Training Data\n');
[acc, p, z] = PredictPheno(BN, false);


AUCWorker(acc, p, z, BN.GetPhenoCol(), true);


