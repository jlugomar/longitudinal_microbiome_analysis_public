function [avgtest, avgdtest, bnf_avgtest_rev, bnf_avgdtest_rev, bnf_avgtest, bnf_avgdtest, data] = CompareDiscretization(loading)
%[testCHAUC, dtestCHAUC, data] = CompareDiscretization()
%
% Compare the accuarcy of a CGBayesNets dataset to BNfinder and Weka
% compare on both mixed Gaussian networks and on networks that have had
% their Gaussian nodes discretized.
%
% for illustration purposes.
%
% 

BNF = true;
bnf_avgtest = 0;
bnf_avgdtest = 0;
bnf_avgtest_rev = 0;
bnf_avgdtest_rev = 0;

generation = 4;
nmaxsims = 9;
n = 20;
ndisc = 5;
pheno = '1';
tsetsize = 25;
verbose = true;

for g = generation:generation + nmaxsims
    fnametrain = ['data_gen', num2str(g),'_train.csv'];
    fnametest = ['data_gen', num2str(g),'_test.csv'];
    if (~loading)
        % generate a network of gaussian variables
        [bn, adjmat] = RandNet(n, ndisc);
        if (sum(adjmat(1,:)) == 0)
            error('no connection to phenotype in random network');
        end
        % simulate that guassian network to generate data
        data = SimNetwork(bn, adjmat);
    else
        % load the dataset:
        data1 = csvread(fnametrain);
        data2 = csvread(fnametest);
        data = [data1;data2];
    end

    names = cell(1,size(data,2));
    disc = IsDiscrete(data);
    for i = 1:size(data,2)
        names{i} = num2str(i);
    end

    % split into training and test:
    traindata = data(1:tsetsize,:);
    testdata = data(tsetsize+1:2*tsetsize,:);
    
    if (~loading)
        csvwrite(fnametrain,traindata);
        csvwrite(fnametest,testdata);
    end
        

    dtraindata = simpleDisc(traindata,10);
    dtestdata = simpleDisc(testdata,10);
    uniquevals = cell(1,n);
    for i = 1:n
        uniquevals{i} = unique([dtraindata(:,i);dtestdata(:,i)],'legacy');
    end
    if (~loading)
        OutputForWEKA(['weka_disc_train', num2str(g),'.arff'], dtraindata, names, uniquevals, true(1,n));
        OutputForWEKA(['weka_disc_test', num2str(g),'.arff'], dtestdata, names, uniquevals, true(1,n));

        OutputForBNFinder(['bnf_weka_train', num2str(g),], traindata, names, disc, false);
        OutputForBNFinder(['bnf_weka_test', num2str(g),], testdata, names, disc, false);

        OutputForBNFinder(['bnf_weka_dtrain', num2str(g),], dtraindata, names, true(1,n), false);
        OutputForBNFinder(['bnf_weka_dtest', num2str(g),], dtestdata, names, true(1,n), false);

        OutputForBNFinder(['bnf_weka_train', num2str(g),], traindata, names, disc, true);
        OutputForBNFinder(['bnf_weka_test', num2str(g),], testdata, names, disc, true);

        OutputForBNFinder(['bnf_weka_dtrain', num2str(g),], dtraindata, names, true(1,n), true);
        OutputForBNFinder(['bnf_weka_dtest', num2str(g),], dtestdata, names, true(1,n), true);
    end

    % BNfinder commands:
    % bnf -e bnf_weka_train2.txt -n out_weka_train2.sif -v -l 3 -c net_weka_train2.cpd
    % bnc -o results_weka_test2.cls -p 1 -c net_weka_train2.cpd -d bnf_weka_test2_test.txt
    %
    % bnf -e bnf_weka_train2_revK2.txt -n out_weka_train2_revK2.sif -v -l 3 -c net_weka_train2_revK2.cpd
    % bnc -o results_weka_test2_revK2.cls -p 1 -c net_weka_train2_revK2.cpd -d bnf_weka_test2_revK2_test.txt
    %
    % bnf -e bnf_weka_dtrain2.txt -n out_weka_dtrain2.sif -v -l 3 -c net_weka_dtrain2.cpd
    % bnc -o results_weka_dtest2.cls -p 1 -c net_weka_dtrain2.cpd -d bnf_weka_dtest2_test.txt
    %
    % bnf -e bnf_weka_dtrain2_revK2.txt -n out_weka_dtrain2_revK2.sif -v -l 3 -c net_weka_dtrain2_revK2.cpd
    % bnc -o results_weka_dtest2_revK2.cls -p 1 -c net_weka_dtrain2_revK2.cpd -d bnf_weka_dtest2_revK2_test.txt
    %
    % check AUCs of BNFinder
    % change to true after running above commands in BNfinder2.0
    % names here with "weka" in them actually refer to BNfinder output
    % files:
    if (BNF)
        wekaresults_gaussian = ['results_weka_test',num2str(g),'_revK2.cls'];
        wekaresults_discretized = ['results_weka_dtest',num2str(g),'_revK2.cls'];
        chauc_gaussian_rev(g) = AUCfromBNfinderFile(wekaresults_gaussian, testdata(:,1));
        chauc_discretized_rev(g) = AUCfromBNfinderFile(wekaresults_discretized, testdata(:,1));
        wekaresults_gaussian = ['results_weka_test',num2str(g),'.cls'];
        wekaresults_discretized = ['results_weka_dtest',num2str(g),'.cls'];
        chauc_gaussian(g) = AUCfromBNfinderFile(wekaresults_gaussian, testdata(:,1));
        chauc_discretized(g) = AUCfromBNfinderFile(wekaresults_discretized, testdata(:,1));
    end
    % then learn a network from that data
    % common parameter values:
    priorPrecision.nu = 1;
    priorPrecision.sigma2 = 1;
    priorPrecision.alpha = 10;
    priorPrecision.maxParents = 3;

    [trainCHAUC(g), ~, testCHAUC(g), ~, ~, ~] = ...
        ModelLearnAndTest(traindata, names, testdata, names, pheno, priorPrecision, ...
            'CompDvG-Gaussian', verbose, names(1:5));

    [dtrainCHAUC(g), ~, dtestCHAUC(g), ~, ~, ~] = ...
        ModelLearnAndTest(dtraindata, names, dtestdata, names, pheno, priorPrecision, ...
            'CompDvG-Discrete', verbose, names);

end

avgtest = sum(testCHAUC(generation:generation+nmaxsims)) / (nmaxsims+1);
avgdtest = sum(dtestCHAUC(generation:generation+nmaxsims)) / (nmaxsims+1);

if (BNF)
    bnf_avgtest_rev = sum(chauc_gaussian_rev(generation:generation+nmaxsims)) / (nmaxsims+1);
    bnf_avgdtest_rev = sum(chauc_discretized_rev(generation:generation+nmaxsims)) / (nmaxsims+1);
    bnf_avgtest = sum(chauc_gaussian(generation:generation+nmaxsims)) / (nmaxsims+1);
    bnf_avgdtest = sum(chauc_discretized(generation:generation+nmaxsims)) / (nmaxsims+1);
end

