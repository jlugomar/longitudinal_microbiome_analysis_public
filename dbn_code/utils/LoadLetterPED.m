function [ids,data] = LoadLetterPED(filename, hasHeaderline)
%[ids,data] = LoadLetterPED(filename, hasHeaderline)
%
% Loads a PED file that has genotypes given by the letters A, C, T, G.
%
% PED format is as follows:
% 1st column is family ID
% 2nd is unique ID.
% 3rd is father ID, 4th is mother ID.
% 5th is gender 1 = male, 2 = female,
% 6th = phenotype; although rarely used.
%
% 7th and 8th are first genotype, etc, going by pairs.
%
% Output datais formatted in a symmetric coding style, where each genotype
% is represented by one number.  
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    hasHeaderline = false;
end

% open file:
fid = fopen(filename);

% first skip the first line if tehre is a header
if (hasHeaderline)
    textscan(fid,'%[^\n]',1);
end

% the pattern given here reads six digits then the rest of the line as a
% string.  It throws away all the digits except the unique ID (col #2)
tcell = textscan(fid,'%*d %d %*d %*d %*d %*d %[\r\n]', 'BufSize', 10000000);
% first cell is the ID field
ids = tcell{1,1};
if (isempty(ids)) % parsing IDs failed because they're string IDs, probably
    % try string identifiers:
    tcell = textscan(fid,'%*s %s %*d %*d %*d %*d %[^\r\n]', 'BufSize', 10000000);
    % first cell is the ID field
    ids = tcell{1,1};
end
fclose(fid);
tlines = tcell{1,2};
clear tcell;

% this will work for both letter codings and numeric codings:
letters = uint8('ACGT1234');

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
data(data > 4) = data(data > 4) -4;

s = size(data);

% each pair of numbers is one alelle.
temphalf = zeros(s(1), s(2)/2,'uint8');
%for i=1:(s(2)/2)
%    temphalf(:,i) = data(:,i*2-1)*5 + data(:,i*2);
%end  

% try doing the above with additions of indexed matrices
odds = (1:(s(2)/2)) * 2 -1;
evens = odds + 1;
temphalf = data(:, odds)*5 + data(:,evens);



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

temphalf(temphalf == 1) = 0;
temphalf(temphalf == 2) = 0;
temphalf(temphalf == 3) = 0;
temphalf(temphalf == 4) = 0;
temphalf(temphalf == 5) = 0;
temphalf(temphalf == 6) = 1;
temphalf(temphalf == 7) = 2;
temphalf(temphalf == 8) = 3;
temphalf(temphalf == 9) = 4;
temphalf(temphalf == 10) = 0;
temphalf(temphalf == 11) = 2;
temphalf(temphalf == 12) = 5;
temphalf(temphalf == 13) = 6;
temphalf(temphalf == 14) = 7;
temphalf(temphalf == 15) = 0;
temphalf(temphalf == 16) = 3;
temphalf(temphalf == 17) = 6;
temphalf(temphalf == 18) = 8;
temphalf(temphalf == 19) = 4;
temphalf(temphalf == 20) = 0;
temphalf(temphalf == 21) = 4;
temphalf(temphalf == 22) = 7;
temphalf(temphalf == 23) = 9;
temphalf(temphalf == 24) = 10;


data = temphalf;


return;
