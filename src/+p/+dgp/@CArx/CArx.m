classdef CArx < handle & p.dgp.CDgp & t.CWithTest;
% autoregressive, i.e.
% x' = \alpha * x + \eps
    % \alpha == 1: unit root
    % \alpha < 1: stationary
    % \alpha > 1: explosive, not implemented
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
    
    properties
        pMu
        pMuLb      % lower bound if random
        pMuUb      % upper bound if random 
        pStd
        pTDistDof  % degree of freedom for the t-dist
        pAlpha
        pDependence
        pHow2Cumulate % for example mAr1
    end
    
    properties(SetAccess = private)
         pVar
    end
    
    methods
        
        %% constructor and setters
        
        function cArx = CArx()
            cArx.pAlpha           = d.CValOptions;
            cArx.pDependence      = d.CValOptions;
            cArx.pMu              = 0;
            cArx.pTDistDof        = 3;
            cArx.pHow2Cumulate    = d.CValOptions('mAr1'); % list additional Dgp Here
            cArx.pHow2Cumulate.pV = 'mAr1'; 
        end
        
        function set.pStd(cArx, v)
            if v > 0
                cArx.pStd = v;
                cArx.pVar = v^2; %#ok<MCSUP>
            else
                error('Std must be positive')
            end
        end
        
        %% simulate data functions
        % done is 2 steps
        % first simulate error term
        % next, cumulate the error term as specified in cArx.pHow2Cumulate
        
        function mSimulate(cArx)
            if isnumeric(cArx.pAlpha.pV) && numel(cArx.pAlpha.pV) == 1 && cArx.pAlpha.pV == 1
                tmp = cArx.pOverSample;
                cArx.pOverSample = 1; % so that the first observation is random
            end
            resDiff  = cArx.mSimulateErrors;
            resLevel = cArx.mAggregate2Level(resDiff);
            cArx.pX  = cArx.mRemoveOverSample(resLevel);
            if isnumeric(cArx.pAlpha.pV) && numel(cArx.pAlpha.pV) == 1 && cArx.pAlpha.pV == 1
                cArx.pOverSample = tmp;
            end
            cArx.mGetFirstDifference();
        end
        
        %% simulate error terms
        
        function res = mSimulateErrors(cArx)
            rw2 = cArx.pT.pV + cArx.pOverSample;
            if cArx.pMu == 0
                indEffect = zeros(rw2, cArx.pN.pV);
            else
                indEffect = a.randBetween(cArx.pMuLb, cArx.pMuUb, rw2, cArx.pN.pV);
            end
            if isempty(cArx.pDependence.pV)
                mu_it = cArx.mIndependentErrors;
            else
                mu_it = cArx.mDependentErrors;
            end
            res = indEffect + mu_it;                   
        end
        
        function eps = mIndependentErrors(cArx)
            rw2 = (cArx.pT.pV + cArx.pOverSample);
            switch cArx.pDistType.pV
                case 't'
                    eps = trnd( cArx.pTDistDof, rw2, cArx.pN.pV);
                otherwise
                    eps  = normrnd(0, cArx.pStd, rw2,  cArx.pN.pV);
            end
        end 
        
        function res = mDependentErrors(cArx)
            switch cArx.pDependence.pV
                case 'wfd'
                    factor_n   = [-0.2, -0.1, -0.3];
                    factor_p   = [0.4,   0.3,  0.2]; 
                    eitTarget  = 0.5;  % how much noise comes from the idiosyncratic term
                case 'sfd'
                    factor_n   = [-0.6, -0.5, -0.75];
                    factor_p   = [0.8,   0.7,  0.7];
                    eitTarget  = 0.2;
                otherwise
                    error('Factor dependence can be wfd or sfd')
            end
            factor_a      = [factor_n, factor_p];
            rw2    = cArx.pT.pV + cArx.pOverSample;
            % error term for the unobserved factors;
            if strcmp(cArx.pDistType.pV, 't')
                theta_t = trnd(cArx.pTDistDof, rw2, numel(factor_n));
            else
                theta_t = normrnd(0, cArx.pStd, rw2,  numel(factor_n));
            end
            % factor loadings
            mu_t = zeros(rw2, cArx.pN.pV);
            for i = 1:cArx.pN.pV
                factInd   = randperm(numel(factor_a), numel(factor_n));  
                f_i       = factor_a(factInd) .* a.randBetween(0.9, 1.1, 1, 3);
                mu_t(:,i) = f_i * theta_t'; 
            end
            % idiosyncratic term
            e_it     = cArx.mIndependentErrors();
            mu_t     =  mu_t - repmat(mean(mu_t), rw2, 1);
            muTarget =  cArx.pStd / mean(std(mu_t));    
            res      =  eitTarget *e_it + (1 -eitTarget) * muTarget * mu_t;
            % disp(cArx.pDistType.pV); disp(mean(res)); disp(std(res));
        end
               
        %% aggregate error terms to get levels
        
        function X = mAggregate2Level(cArx, eps)
            if numel(cArx.pAlpha.pV) == 1
                X = cArx.(cArx.pHow2Cumulate.pV)(eps);
            else % mixed
                pAlphaM  = cArx.pAlpha.pV;
                n        = cArx.pN.pV / numel(pAlphaM);
                if n ~= floor(n)
                    error('N must be divisible with the number of alphas.')
                end
                rw2 = cArx.pT.pV + cArx.pOverSample;
                X = zeros(rw2, cArx.pN.pV);
                tmpN           = cArx.pN.pV; 
                cArx.pN.pV     = n;
                for i = 1:numel(pAlphaM)
                    if iscell(pAlphaM(i))
                        cArx.pAlpha.pV = pAlphaM{i};
                    else
                        cArx.pAlpha.pV = pAlphaM(i);
                    end
                    cl_0 = (i-1)* n +1;
                    cl_1 = i    * n;
                    X(:, cl_0: cl_1) = cArx.(cArx.pHow2Cumulate.pV)(eps(:, cl_0: cl_1));
                end
                cArx.pAlpha.pV = pAlphaM;
                cArx.pN.pV     = tmpN; 
            end
        end
        
        function X = mAr1(cArx, eps)
            % combine
            X      = zeros(size(eps));
            X(1,:) = eps(1,:);
            for t = 2:u.rw(eps)
                X(t, :) = cArx.pAlpha.pV .* X(t-1,:) + eps(t, :);
            end
        end
        
        function res = mRemoveOverSample(cArx, X)
            tmp     = X((cArx.pOverSample+1):end,:);                % throw away first few obserations
            tmpInit = repmat(X(cArx.pOverSample+1,:),u.rw(tmp),1);  % start at zero          
            res     = tmp - tmpInit;
        end
        
        %% iterate functions
        % each iterate function serves to go over 1 or more combinations
        % once it is finished, it can call additional iterate functions
        % eventuelly the chain has to reach the mDoSimulation function
        % mDoSimulation will call obj.(methodStr) obj.cArx.pSimCount times
        % cl should be incremented by the number of columns the given iterate occupies
        % for i = 1:par1
        %    for j = 1:par2
        %        cArx.mIterate(obj, methodStr, cl+2)
        
        function mIterateWithDistDependence(cArx, obj, methodStr)
            cArx.mInitIterations();
            for i = 1:numel(cArx.pDependence.pOptions)
                cArx.mAddParameter('pDependence', i) 
                cArx.mIterateWithDist(obj, methodStr);
            end
        end
        
        function mIterateWithDist(cArx, obj, methodStr)
            cArx.mInitIterations();
            for i = 1:numel(cArx.pDistType.pOptions)
                cArx.mAddParameter('pDistType', i) 
                cArx.mIterate(obj, methodStr);
            end
        end
            
        function mIterateNoDist(cArx, obj, methodStr)
            cArx.mInitIterations();
            cArx.mIterate(obj, methodStr);
        end
        
        function res = mIterateTest(cArx, dgp, res, rw) %#ok<INUSL>
            res(rw, 1) = normrnd(0, 1);
            res(rw, 2) = normrnd(1, 1);
        end
        
    end
        
    methods(Access = protected)      
        
        function mIterate(cArx, obj, methodStr)
            for i = 1:numel(cArx.pAlpha.pOptions)
                cArx.mAddParameter('pAlpha', i);
                for j = 1:numel(cArx.pT.pOptions)
                    cArx.mAddParameter('pT', j);
                    for k = 1:numel(cArx.pN.pOptions)
                        cArx.mAddParameter('pN', k);
                        cArx.mDoSimulations(obj, methodStr);
                    end
                end    
            end
        end
       
    end
    
    methods(Static)
        
        %% test DGP
        % create an instance
        % set parameters
        % simulate
        % plot
        
        function cArx = mCreateObj()
            cArx = p.dgp.CArx;
            cArx.pHow2Cumulate.pV = 'mAr1';
            cArx.pStd             = 1;
            cArx.pSimCount        = 100;
            cArx.pOverSample      = 200;
            cArx.pT.pOptions      = [100, 200, 500];
            cArx.pT.pV            = 100;
            cArx.pN.pOptions      = [1, 4];
            cArx.pN.pV            = 4;
            cArx.pAlpha.pOptions  = [1, 0.9];
            strPathTmp            = mfilename('fullpath');
            [strPath, ~, ~]       = fileparts(strPathTmp);
            cArx.pExportDir       = [strPath, filesep, 'testRes'];
        end
        
        function [cArx] = mTestUnitRootIndependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = 1;
            cArx2.pAlpha.pV = 1;
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestUnitRootWeakDependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = 1;
            cArx.pDependence.pV  = 'wfd';
            cArx2.pAlpha.pV = 1;
            cArx2.pDependence.pV = 'wfd';
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestUnitRootStrongDependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = 1;
            cArx.pDependence.pV  = 'sfd';
            cArx2.pAlpha.pV = 1;
            cArx2.pDependence.pV = 'sfd';
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestStationaryIndependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = 0.9;
            cArx2.pAlpha.pV = 0.9;
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestStationaryWeakDependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = 0.9;
            cArx.pDependence.pV  = 'wfd';
            cArx2.pAlpha.pV = 0.9;
            cArx2.pDependence.pV = 'wfd';
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestStationaryStrongDependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = 0.9;
            cArx.pDependence.pV  = 'sfd';
            cArx2.pAlpha.pV = 0.9;
            cArx2.pDependence.pV = 'sfd';
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
               
        function [cArx] = mTestMixedIndependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = [1, 0.9];
            cArx2.pAlpha.pV = [1, 0.9];
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestMixedWeakDependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = [1, 0.9];
            cArx.pDependence.pV  = 'wfd';
            cArx2.pAlpha.pV = [1, 0.9];
            cArx2.pDependence.pV = 'wfd';
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx] = mTestMixedStrongDependent()
            [cArx, cArx2]   = p.dgp.CArx.mCreateTwoCArx();
            cArx.pAlpha.pV  = [1, 0.9];
            cArx.pDependence.pV  = 'sfd';
            cArx2.pAlpha.pV = [1, 0.9];
            cArx2.pDependence.pV = 'sfd';
            p.dgp.CArx.mSimulatePlotResult(cArx, cArx2)
            res = true;
        end
        
        function [cArx, cArx2] = mCreateTwoCArx()
            cArx = p.dgp.CArx.mCreateObj();
            cArx.pDistType.pV  = 'n';
            cArx2 = p.dgp.CArx.mCreateObj();
            cArx2.pDistType.pV = 't';
        end
        
        function mSimulatePlotResult(cArx, cArx2)
            cArx.mSimulate();
            cArx2.mSimulate();
            cc = floor(a.randBetween(1, u.cl(cArx.pX), 1, 1));
               
            subplot(2,1,1);
            plot(cArx.pX(:,cc));
            hold on
            plot(cArx2.pX(:,cc));
            title('Levels');
            legend('normal','student-t');
            
            subplot(2,1,2);
            plot(cArx.pDX(:,cc));
            hold on
            plot(cArx2.pDX(:,cc));
            title('First differences');
            legend('normal','student-t');
            
            corrplot(cArx.pDX );
            title('Normal distribution');
            corrplot(cArx2.pDX);
            title('Student-t distribution');
            
        end
        
        %% test iterate functions
        
        function [cArx] = mTestIterateNoDist()
            cArx = p.dgp.CArx.mCreateObj();
            cArx.pFileName = ['mTestIterateNoDist', '.csv'];
            cArx.pHeaderOutput = {'a', 'b'};
            cArx.mIterateNoDist(cArx, 'mIterateTest')
            cArx.mExportRes();
        end

        function [cArx] = mTestIterateWithDist()
            cArx = p.dgp.CArx.mCreateObj();
            cArx.pDistType.pOptions = {'n', 't'};
            cArx.pAlpha.pOptions = {1, 0.9, {1, 0.9}};
            cArx.pN.pOptions     = [2, 4];
            cArx.pHeaderOutput = {'c', 'd'};
            cArx.mIterateWithDist(cArx, 'mIterateTest')
            cArx.pFileName = ['mTestIterateWithDist', '.csv'];
            cArx.mExportRes();
        end
   
        function [cArx] = mTestIterateWithDistDependence()
            cArx = p.dgp.CArx.mCreateObj();
            cArx.pFileName = ['mTestIterateWithDistDependence', '.csv'];
            cArx.pDistType.pOptions = {'t', 'n'};
            cArx.pDependence.pOptions = {'wfd', 'sfd'};
            cArx.pHeaderOutput = {'e', 'f'};
            cArx.mIterateWithDistDependence(cArx, 'mIterateTest')
            cArx.mExportRes();
        end
                   
    end
    
end