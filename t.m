% this project tests classes in src
%% init
clear
close all
rng('default')
tic;

%% data generating process
cDgp   = p.dgp.CArx;
cDgp.pExportDir = 'C:\Users\Leni\Dropbox\Sajat\Research\03_Projects\2015_Science_BCCStatisticalTesting_UR_FPML\4_src_2017\test';
cDgp.mRunMethodsStarting('mTest');

%% boudary crossing representation
cCount = p.bcc.CCount;
cDgp.mRunMethodsStarting('mTest');

%% unit root tests
cUnitRootTests = p.s.CUnitRootTests;
cDgp.mRunMethodsStarting('mTest');

%% measure runTime
toc;