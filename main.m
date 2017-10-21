%% init
clear
close all
rng('default')
warning('off');
strPathTmp = mfilename('fullpath');
[strPath, ~, ~] = fileparts(strPathTmp);
if isempty(strPath)
    warning('mfilename only works if the script is called from command line. For debugging, pls. fill out strPath manually.')
end
addpath([strPath, filesep, 'src'])
tic;

%% 
cDgp            = p.dgp.CArx;
cDgp.pExportDir = [strPath, filesep,'res'];   
isTest = false;
if isTest
    cDgp.pSimCount = 100;
else
    cDgp.pSimCount = 2000; %#ok<UNRCH>
end
cDgp.pStd           = 1;
cDgp.pOverSample    = 200;

%% how boundaries effect bcc power
cDgp.pAlpha.pOptions = [1, 0.9];
cDgp.pN.pOptions     = 1;
cDgp.pT.pOptions     = [100, 500, 1000];
cBcc               = p.bcc.CCount;
cBcc.pL.pOptions   = (-1:-1:-6);
cBcc.pU.pOptions   = - cBcc.pL.pOptions;
cBcc.pRestart      = 0;
cDgp.pHeaderOutput = {'b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6'};
cDgp.mIterateNoDist(cBcc, 'mIterateOverBoundaries');
cDgp.pFileName = 'effectOfBoundaries.csv';
cDgp.mExportRes;

%% comparative mc unit root tests
cUnitRoot           = p.s.CUnitRootTests();

%% for time series
cDgp.pDistType.pOptions = {'n', 't'};
cDgp.pAlpha.pOptions    = [1, 0.9];
cDgp.pN.pOptions        = 1;
cDgp.pT.pOptions        = [50,100,200];
cDgp.pHeaderOutput      = {'ADF', 'PP', 'VR', 'BCC'};
cDgp.mIterateWithDist(cUnitRoot, 'mTimeSeriesTests');
cDgp.pFileName = 'timeSeriesUnitRootTests.csv';
cDgp.mExportRes;

%% for first generation panel unit root
cDgp.pDistType.pOptions = {'n', 't'};
cDgp.pAlpha.pOptions    = {1, 0.9, {1, 0.9}};
cDgp.pN.pOptions        = [12, 20];
cDgp.pT.pOptions        = [25, 50];
cDgp.pHeaderOutput      = {'IPS_1', 'IPS_2', 'IPS_3', 'MW_1', 'MW_2' 'CH_1', 'CH_2', 'BCC'};
cDgp.mIterateWithDist(cUnitRoot, 'mPanelFirstGenTests');
cDgp.pFileName = 'firstGenPanelUnitRootTests.csv';
cDgp.mExportRes;

%% measure runTime
toc;