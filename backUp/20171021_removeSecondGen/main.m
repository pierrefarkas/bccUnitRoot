%% init
clear
close all
rng('default')
strPath = mfilename('fullpath');
addpath([strPath, filesep, 'src'])
tic;

cDgp            = p.s.dgp.CArx;
cDgp.pExportDir = [strPath, filesep,'res'];   
isTest = true;
if isTest
    cDgp.pSimCount = 100;
else
    cDgp.pSimCount = 2000; %#ok<UNRCH>
end

cDgp.pStd           = 1;
cDgp.pOverSample    = 200;
cDgp.pT.pOptions     = [100, 500, 1000];
cDgp.pN.pOptions     = 1;
cDgp.pAlpha.pOptions = [1, 0.9];

%% how boundaries effect bcc power
cBcc               = p.bcc.CCount;
cBcc.pL.pOptions   = (-1:-1:-6);
cBcc.pU.pOptions   = - cBcc.pL.pOptions;
cBcc.pRestart      = 0;
cDgp.pHeaderOutput = {'[-1, 1]', '[-2, 2]', '[-3, 3]', '[-4, 4]', '[-5, 5]', '[-6, 6]'};
cDgp.mIterateNoDist(cBcc, 'mIterateBoundaries');
cDgp.mExportRes;

%% comparative mc unit root tests
cUnitRoot           = p.s.CUnitRootTests();

%% for time series
cDgp.pT.pOptions   = [50,100,200];
cDgp.pHeaderOutput = {'ADF', 'PP', 'VR', 'BCC'};
cDgp.mIterateWithDist(cUnitRoot, 'mTimeSeriesTests');
cDgp.mExportRes;

%% for first generation panel unit root
cDgp.pT.pOptions     = [25, 50];
cDgp.pN.pOptions     = [12, 20];
cDgp.pAlpha.pOptions = [1, 0.9, [1, 0.9]];
mcPanelFirstGen     = zeros(cDgp.pCaseCount, (cDgp.pClCount + cUnitRoot.pFirstGenPanelCount));
headers             = {'Dist', '\alpha', 'N', 'T', 'IPS-1', 'IPS-2', 'IPS-3', 'MW-1', 'MW-2' 'CH-1', 'CH-2', 'BCC'};
cDgp.mIterateWithDist(mcPanelFirstGen, cUnitRoot, 'mPanelFirstGenTests');
u.csvwrite_with_headers('.\res\mcPanelFirstGen.csv', mcPanelFirstGen, headers);

%% second generation panel unit root
cDgp.pAlpha.pOptions = ['wfd', 'sfd'];  % weak factor dependence, strong factor dependence
mcPanelSecondGen    = zeros(cDgp.pCaseCount, (cDgp.pClCount + cUnitRoot.pSecondGenPanelCount));
headers             = {'Dist', '\alpha', 'N', 'T', 'BNG-1', 'BNG-2', 'PS-1', 'PS-2', 'BCC'};
cDgp.mIterateWithDistDependence(mcPanelSecondGen, cUnitRoot, 'mPanelSecondGenTests');
u.csvwrite_with_headers('.\res\mcPanelSecondGen.csv', mcPanelSecondGen, headers);