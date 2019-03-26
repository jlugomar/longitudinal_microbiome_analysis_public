%% bnpathscript
% scriptfile for setting up CGBayesNets path variables
%
% must be run from the CGBayesNets home directory 
addpath(pwd);
addpath([pwd,filesep,'data']);
addpath([pwd,filesep,'data',filesep,'internal']);
addpath([pwd,filesep,'data',filesep,'testdata']);
addpath([pwd,filesep,'data',filesep,'testdata',filesep,'internal']);
addpath([pwd,filesep,'discretization']);
addpath([pwd,filesep,'utils']);
addpath([pwd,filesep,'utils',filesep,'internal']);
addpath([pwd,filesep,'auctools']);
addpath([pwd,filesep,'auctools',filesep,'internal']);
addpath([pwd,filesep,'inference']);
addpath([pwd,filesep,'inference',filesep,'internal']);
addpath([pwd,filesep,'netlearning']);
addpath([pwd,filesep,'netlearning',filesep,'internal']);

warning('off','all')
