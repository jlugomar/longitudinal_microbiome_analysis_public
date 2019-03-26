function ProcessGeoSet(norm)
% load the GSE19301 dataset and process the data.
%
% this does the following:
% 1) remove unnecessary data columns
% 2) use AffyGeneMap_GPL96-15653_simple.txt to rename mRNA probes with the
% cannonical gene names.
% 3) average the non-exacerbation samples per individual, and then subtract
% this from the one exacerbation sample to get a differential expression
% measurement.
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

if (nargin < 1)
    norm = true;
end

[phndata, phncols, ~, ~, ~] = RCSVLoad('GSE19301_pheno_matrix_simple.txt',false,'\t');
[numdata, cols, strdata, ~, ~] = RCSVLoad('GSE19301_data_matrix.txt',false,'\t');
genenames = strdata{1}';
genedata = numdata';
sample_ids = cols(2:end)';

% also log normalize each gene expression column:
if (norm)
    for i = 1:length(genenames)
        genedata(:,i) = (log(genedata(:,i)) - mean(log(genedata(:,i)))) ./ std(log(genedata(:,i)));
    end
end


pheno = 'Exacerbation';

dropcols = {'any relevant respiratory infections','atopy status',...
    'any inhaled cs use  0=no  1=yes','any intranasal cs use  0=no  1=yes',...
    'any systemic cs use  0=no  1=yes','visit number', 'ID'};

drops = zeros(size(phncols));
for i = 1:length(dropcols)
    drops = drops | strcmpi(dropcols(i),phncols);
end

% drop the columns we're not interested in
phncols = phncols(~drops);
phndata = phndata(:,~drops);

% move phenotype to column #1
phenocol = strcmpi(pheno,phncols);
phndata = [phndata(:,phenocol),phndata(:,~phenocol)];
phncols = {pheno, phncols{~phenocol}};

% replace Affy ID's with actual gene symbols, where available -
[~, genecols, genestrdata, ~, ~] = RCSVLoad('AffyGeneMap_GPL96-15653_simple.txt',false,'\t');
% probes with no gene id listed
probe_ids = genestrdata{1,1};
gene_probenames = genestrdata{1,3};
UNKNOWN = 'unknown'; 
genenamemap = containers.Map('KeyType','char','ValueType','int32');
for i = 1:length(genenames)
    % see if this probe is in the probe-gene map
    matchinds = strcmpi(genenames(i),probe_ids);
    % if the name is "unknown" just use the affy probe name
    if ((sum(matchinds) > 0) && (sum(strcmpi(UNKNOWN,gene_probenames{matchinds})) == 0))
        firstmatch = find(matchinds);
        probename = gene_probenames(firstmatch(1));
        probename = probename{:};
        if (isKey(genenamemap,probename))
            % we've seen this gene name before, so increment number of
            % probes seen:
            numseen = genenamemap(probename);
            numseen = numseen + 1;
        else
            % otherwise, first time we've seen this gene name
            numseen = 1;            
        end
        genenames{i} = [probename,'_p',num2str(numseen)];
        % update quantity
        genenamemap(probename) = numseen;
    end
end


allcols = {phncols{:}, genenames{:}};
alldata = [phndata, genedata];



%% make this into a differential expression dataset:
% find all the expressions for each person:
DONOR = 'Donor';
donorcol = find(strcmpi(DONOR,phncols)); 
% since the gene cols are stacked after the phenotypic/demographic cols, 
% the index is the same for DONORCOL 
donors =  unique(phndata(:,donorcol));

% don't subtract out the race column:
RACE = 'race';
racecol = find(strcmpi(RACE,phncols));

diffs = zeros(length(donors),length(allcols));
for i = 1:length(donors)
    dset = phndata(:,donorcol) == donors(i);
    exset = dset & phndata(:,1) == 1;
    if (sum(exset) > 0)
        baseset = dset & phndata(:,1) == 0;
        basemeans = mean(alldata(baseset,:));
        if (sum(exset) > 1)
            diffs(i,:) = mean(alldata(exset,:)) - basemeans;
        else
            diffs(i,:) = alldata(exset,:) - basemeans;
        end
        diffs(i,racecol) = mean(alldata(exset,racecol));
    end
end

% save file:
TitlesAndTabsFileout('GSE19301_diffexpr_geneids.txt', allcols, diffs);






