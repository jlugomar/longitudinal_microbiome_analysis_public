function [ids,snpdata,rsname,chr,position,majoralleles] = LoadGWAS(pedfile, mapfile, hasHeaderline, outputSymmetric, cleaning,forceStringID)
% [ids,snpdata,rsname,chr,position,majoralleles] = LoadGWAS(pedfile, mapfile, hasHeaderline, outputSymmetric, cleaning,forceStringID)
%
% function to load SNP data from PED/MAP pair of files.
% load SNP names from MAP file
% uses efficient coding and loading to read a PED file of (1,2,3,4)-coded
% genotypes into a matrix of UINT8's.  Uses (0,1,2) minor allele-count
% coding.
%
% updated: filters out all 3rd, 4th etc. most popular genotypes now; which
% makes it consistent with coding output by plink --recode 12, where 1 is
% the minor allele, 2 is the major allele, and 0 is a missing genotype.
%
% Columns 7th and 8th are first genotype, etc, going by pairs.
%
% INPUT:
%   filename : a plink-formatted PED file
%   mapfile : opional, a plink formatted MAP file that describes variants
%       in the PED file
%   hasHeaderline : optional, default = false.
%   outputSymmetric : if true, will output in 1-10 coding for combinations
%       of specific alleles; if false, will output in (0,1,2) minor allele count
%       coding. default = false (0,1,2)-coding.
%   cleaning : boolean (optional).  If true, will remove SNPs with fewer
%       than 10% coded values
%
% OUTPUT:
%   ids : the subject IDs; parallel array to the snpdata
%   snpdata : actual SNP coded values for each subject x SNP
%   rsname : list of SNP names, parallel to SNPDATA
%   chr : chromosome nubmers of each SNP, parallel array to RSNAME
%   position : SNP chromosomal positions, parallel array to RSNAME
%   majorallels : list of major allels for each SNP
%
% 
% Copyright Michael McGeachie, 2014.  MIT license. See cgbayesnets_license.txt.

skipmap = false;
if (nargin < 6)
    forceStringID = false;
end
if (nargin < 5)
    cleaning = true;
end
if (nargin < 4)
    outputSymmetric = false;
end
if (nargin < 3)
    hasHeaderline = false;
end
if (nargin < 2 || isempty(mapfile))
    skipmap = true;
    rsname = {};
    chr = [];
    position = [];
end

% load SNP names from MAP file
if (~skipmap)
    fid = fopen(mapfile);
    if (hasHeaderline)
        textscan(fid,'%[^\r\n]',1);
    end

    % this picks out the chromosome number and the rs name (snp name)
    % and the chromosome position
    mcell = textscan(fid,'%d %s %*d %d', 'BufSize', 10000000);
    fclose(fid);
    chr = mcell{1};
    rsname = mcell{2};
    position = mcell{3};
    clear mcell;
end

% open PED file:
fid = fopen(pedfile);
% first skip the first line if there is a header
if (hasHeaderline)
    textscan(fid,'%[^\r\n]',1);
end

% the pattern given here reads six digits then the rest of the line as a
% string.  It throws away all the digits except the unique ID (col #2)
tcell = textscan(fid,'%*d %d %*d %*d %*d %*d %[^\r\n]', 'BufSize', 10000000);
% first cell is the ID field
ids = tcell{1,1};
if (forceStringID || isempty(ids)) % parsing IDs failed because they're string IDs, probably
    % try string identifiers:
    fclose(fid);
    fid = fopen(pedfile);
    if (hasHeaderline)
        textscan(fid,'%[^\r\n]',1);
    end        
    tcell = textscan(fid,'%*s %s %*d %*d %*d %*d %[^\r\n]', 'BufSize', 10000000);
    % first cell is the ID field
    ids = tcell{1,1};
end
fclose(fid);
tlines = tcell{1,2};
clear tcell;

% this will work for both letter codings and numeric codings:
letters = uint8('ACGT01234');

% pre-allocate data:
dline = uint8(tlines{1});
% tricky trick here is actually casting every character of the PED file
% into a uint8, which includes the spaces.  So we skip those with this
% staggered index:
inds = [1:2:length(dline)];
data = zeros(length(ids),length(inds),'uint8');

% cast each string to uints and stack it up together:
for i = 1:length(tlines)
    dline = uint8(tlines{i});
    data(i,:) = dline(inds);
end

clear tlines dline inds;

% now need to parse each letter and turn it into a number.
for l = 1:length(letters)
    data(data == letters(l)) = l;
end
data(data > uint8(4)) = data(data > uint8(4)) - uint8(5);

s = size(data);
odds = (1:(s(2)/2)) * 2 -1;
evens = odds + 1;

if (~outputSymmetric)
    % count occurrences of each allele:
    letters = uint8([0,1,2,3,4]);
    genotypes = zeros(length(letters),s(2));
    for l = 1:length(letters)
        genotypes(l,:) = sum(data == letters(l));
    end
    % use additions of indexed matrices
    genocounts = genotypes(:,odds) + genotypes(:,evens);
    [~,majoralleles] = max(genocounts);
    singleinds = 1:length(majoralleles);
    singleinds = (singleinds -1) * size(genotypes,1) + majoralleles;
    genocounts(singleinds) = 0;
    [~,minoralleles] = max(genocounts);
    minoralleles = letters(minoralleles);
    majoralleles = letters(majoralleles);
    if (cleaning)
        skip = genocounts(1,:) > (0.1 * sum(genocounts));
        if (~isempty(rsname))
            rsname = rsname(~skip);
        end
        if (~isempty(chr))
            chr = chr(~skip);
        end
        if (~isempty(position))
            position = position(~skip);
        end
    end
    minoralleles = uint8(minoralleles);
    snpdata = zeros(s(1),s(2)/2,'uint8');
    for i=1:s(2)/2
        % oddly, this does not seem to be very slow:
        snpdata(:,i) = (data(:, odds(i)) == minoralleles(i)) + (data(:,evens(i)) == minoralleles(i));
    end
    if (cleaning)
        snpdata = snpdata(:,~skip);
        majoralleles = majoralleles(:,~skip);
    end

else
    majoralleles = [];
    snpdata = data(:, odds)*5 + data(:,evens);

    % symmetric coding
    % 1,1 = 1 
    % 1,2 = 2
    % 1,3 = 3 
    % 1,4 = 4
    % 2,1 = 2
    % 2,2 = 5
    % 2,3 = 6
    % 2,4 = 7
    % 3,1 = 3
    % 3,2 = 6
    % 3,3 = 8
    % 3,4 = 9
    % 4,1 = 4
    % 4,2 = 7
    % 4,3 = 9
    % 4,4 = 10

    snpdata(snpdata == 1) = 0;
    snpdata(snpdata == 2) = 0;
    snpdata(snpdata == 3) = 0;
    snpdata(snpdata == 4) = 0;
    snpdata(snpdata == 5) = 0;
    snpdata(snpdata == 6) = 1;
    snpdata(snpdata == 7) = 2;
    snpdata(snpdata == 8) = 3;
    snpdata(snpdata == 9) = 4;
    snpdata(snpdata == 10) = 0;
    snpdata(snpdata == 11) = 2;
    snpdata(snpdata == 12) = 5;
    snpdata(snpdata == 13) = 6;
    snpdata(snpdata == 14) = 7;
    snpdata(snpdata == 15) = 0;
    snpdata(snpdata == 16) = 3;
    snpdata(snpdata == 17) = 6;
    snpdata(snpdata == 18) = 8;
    snpdata(snpdata == 19) = 4;
    snpdata(snpdata == 20) = 0;
    snpdata(snpdata == 21) = 4;
    snpdata(snpdata == 22) = 7;
    snpdata(snpdata == 23) = 9;
    snpdata(snpdata == 24) = 10;
    
end


return;
