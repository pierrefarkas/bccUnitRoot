classdef CUnitRootTests < handle & t.CWithTest; 
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com 
    
    properties(Constant)
        pTsCount             = 4;
        pFirstGenPanelCount  = 8;
        pSecondGenPanelCount = 5;
    end
    
    properties
        pAlpha % by default 5%
    end
    
    methods
        
        function cTest = CUnitRootTests()
            cTest.pAlpha = 0.05;
        end
        
        function res = mTimeSeriesTests(cTest, cDgp, res, rw)
            res(rw, 1) = adftest   (cDgp.pX, 'alpha', cTest.pAlpha, 'model', 'ARD' ) + 0;    % Augmented Dickey Fuller, individual effect under the alternative is non-zero 
            res(rw, 2) = pptest    (cDgp.pX, 'alpha', cTest.pAlpha, 'model', 'ARD' ) + 0;    % Phillip Peron Test, individual effect under the alternative is non-zero
            % res(rw, 3) = kpsstest   (cDgp.pX, 'alpha', 0.05, 'trend', false ) + 0;    % KPSS Test
            res(rw,3)  = vratiotest(cDgp.pX, 'alpha', cTest.pAlpha, 'IID', true ) + 0;       % Variance Ratio Test
            cBcc          = p.bcc.CCount();
            cBcc.pDX      = cDgp.pDX;
            cBcc.mUnitRootFirstGenOptimalBoundary();
            res(rw, 4) = cBcc.pD;
        end
        
        function res = mPanelFirstGenTests(cTest, cDgp, res, rw)
           % try
                IPS        = p.s.unitRoot.Test_IPS(cDgp.pX, 1);
                res(rw, 1) = IPS.Wbar_pvalue    < cTest.pAlpha;
                res(rw, 2) = IPS.Zbar_pvalue    < cTest.pAlpha;
                res(rw, 3) = IPS.Zbar_DF_pvalue < cTest.pAlpha;
            % catch
            %    res(rw, 1) = NaN;
            %    res(rw, 2) = NaN;
            %    res(rw, 3) = NaN;
            % end
            % try
                MW1 = p.s.unitRoot.Test_MW(cDgp.pX, 2);         
                res(rw, 4) = MW1.PMW_pvalue < cTest.pAlpha;   % Pvalue Pooled test statistic (Maddala Wu 1999)
                res(rw, 6) = MW1.ZMW_pvalue < cTest.pAlpha;   % Pvalue Pooled test statistic (Choi 2001)
            % catch
            %    res(rw, 4) = NaN;
            %    res(rw, 6) = NaN;
            % end
            % try
                MW2 = p.s.unitRoot.Test_MW(cDgp.pX, 2, NaN,0);         
                res(rw, 5) = MW2.PMW_pvalue < cTest.pAlpha;   % Pvalue Pooled test statistic (Maddala Wu 1999)
                res(rw, 7) = MW2.ZMW_pvalue < cTest.pAlpha;   % Pvalue Pooled test statistic (Choi 2001)
            % catch
            %    res(rw, 5) = NaN;
            %    res(rw, 7) = NaN;
            % end
            cBcc          = p.bcc.CCount();
            cBcc.pDX      = cDgp.pDX;
            cBcc.mUnitRootFirstGenOptimalBoundary();
            res(rw, 8)    = cBcc.pD;
        end
        
        % not yet implemented
        % function res = mPanelSecondGenTests(cTest, cDgp, res, rw)
            % Choi(2002 ) Test of Unit Root "Combinaition Unit Root Tests for Cross-Sectionally Correlated Panels"
        %    try
        %        CH  = Test_ch_f(dgp.data,2,zeros(1,dgp.n(c_n)),0);
        %        res(1,c_nsim) = CH.Pm_pvalue;             % Pvalue associated to the Fisher's modified Statistic
        %        res(2,c_nsim) = CH.Z_pvalue;              % Pvalue associated to the Inverse Normal Test
        %        res(3,c_nsim) = CH.Lstar_pvalue;          % Pvalue associated to the Modified Logit Test
        %    catch
        %        res(1,c_nsim) = NaN;                      % Pvalue associated to the Fisher's modified Statistic
        %        res(2,c_nsim) = NaN;                      % Pvalue associated to the Inverse Normal Test
        %        res(3,c_nsim) = NaN;                      % Pvalue associated to the Modified Logit Test
        %    end
            % Pesaran (2004) unit root test "A simple panel unit root test in the presence of cross section dependence", Mimeo
        %    try
        %        PS = Test_ps_f(dgp.data,2,0,0);
        %        res(4,c_nsim) = PS.CIPS_pvalue;         % Pvalue of the CIPS statistic
        %        res(5,c_nsim) = PS.CIPS_star_pvalue;    % Pvalue of the truncated CIPS* statistic
        %    catch
        %        res(4,c_nsim) = NaN;        % Pvalue of the CIPS statistic
        %        res(5,c_nsim) = NaN;        % Pvalue of the truncated CIPS* statistic
        %    end
            % Not yet implemented
            % cBcc  =  p.bcc.CCount();
            % cBcc
            % cBcc.mUnitRootDependent(cDgp.pDx);
        %    res(rw, 6) = NaN;   
        % end
        
    end
  
    methods(Static, Access = private)
    
        function [cDgp] = mTestUnitRoot()            
            % 1. dgp 2. analytics 3. simulate
            cDgp                = p.dgp.CArx;
            cDgp.pExportDir     = 'C:\Users\Leni\Dropbox\Sajat\IT\Development\Matlab\Active\+p\+s\@CUnitRootTests\res'; 
            cDgp.pSimCount      = 50;
            cDgp.pStd           = 1;
            cDgp.pOverSample    = 200;
            % list options in order
            cDgp.pDistType.pOptions = {'n', 't'};
            cDgp.pAlpha.pOptions    = [1, 0.9];
            % 2. analytics
        end
    
    end
    
    methods(Static)
       
        function [cDgp] = mTestTimeSeries()
            cDgp = p.s.CUnitRootTests.mTestUnitRoot();
            cDgp.pN.pOptions   = 1;
            cDgp.pT.pOptions   = [50, 100, 200];
            cDgp.pHeaderOutput = {'ADF', 'PP', 'VR', 'BCC'};
            cUnitRoot          = p.s.CUnitRootTests();
            cDgp.mIterateWithDist(cUnitRoot, 'mTimeSeriesTests');
            cDgp.mExportRes();
        end
       
        function [cDgp] = mTestFirstGenPanel()
            cDgp = p.s.CUnitRootTests.mTestUnitRoot();
            cDgp.pN.pOptions   = [12, 20];
            cDgp.pT.pOptions   = [25, 50];
            cDgp.pHeaderOutput = {'IPS_wBar', 'IPS_tBar', 'IPS_zBar', 'MW_lag0', 'MW_lagDF', 'CH_lag0', 'CH_lagDF', 'BCC'};
            cUnitRoot          = p.s.CUnitRootTests();
            cDgp.mIterateWithDist(cUnitRoot, 'mPanelFirstGenTests');
            cDgp.mExportRes();
        end
        
    end
    
end