classdef CCount < handle & t.CWithTest; 
% Non-parametric boundary crossing counting
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
    
    % input data
    properties
        pDX
        pStd
        rw
        cl
        pX     % to explain pDX
    end
    
    % count settings
    properties
       pL       % lower bound
       pU       % upper bound
       pRestart % how to restart, 0 for zero, empty otherwise
    end
    
    % counting outcome
    properties
        pTuc   % total upper crossing
        pTlc   % total lower crossing
        pTbc   % total boundary crossing (upper and lower both)
        pP     % upper crossing probability
        pQ     % lower crossing probability
        pFET   % first exit time 
        pBCA   % boundary crossing events
        pBCD   % upper minus lower crossing
        pLcc   % lower crossing count
        pUcc   % upper crossing count
        cTable % contingency table
        pPc    % convergence probability 
    end
    
    % unit test outcome
    properties
        pPVal    % p-value
        pAlpha   % significance level
        pD       % accept: 0, reject: 1   
        
    end
    
    methods
        
        function cCount = CCount(dX)
            if nargin > 0
                cCount.pDX = dX;
            end
            cCount.pL       = d.CValOptions;
            cCount.pU       = d.CValOptions;
            cCount.pRestart = 0;
            cCount.pAlpha   = 0.05;
        end
        
        function set.pDX(cCount, dX)
            cCount.pDX  = dX;
            cCount.mSetDx()
        end
        
        %% counting
        % set boundary or use optimal bounday
        % represent the data using bcc events
        % classify the events based on pX
        
        function mOptimalBoundary(cCount)
            if cCount.rw < 100 
                cCount.pU.pV = 1;
            else
                cCount.pU.pV = 1 + min(1, 1 / cCount.cl) * ((cCount.rw - 100)/225);
            end
            cCount.pL.pV = - cCount.pU.pV; 
        end
        
        function mCount(cCount)            
            % init vars
            lr          = cCount.pDX; 
            L           = cCount.pL.pV * cCount.pStd;
            U           = cCount.pU.pV * cCount.pStd;
            cCount.pLcc = zeros(cCount.rw+1, cCount.cl);
            cCount.pUcc = zeros(cCount.rw+1, cCount.cl);
            cCount.pBCA = zeros(cCount.rw+1, cCount.cl);
            cCount.pBCD = zeros(cCount.rw+1, cCount.cl);
            cCount.pFET = zeros(cCount.rw+1, cCount.cl); 
            lra         = zeros(cCount.rw+1, cCount.cl);  % restarted
            t           = zeros(cCount.rw+1, cCount.cl);  % time
            zerov       = zeros(1,cCount.cl);             % zero vector  
            onev        = ones(1,cCount.cl);              % one vector
            switch cCount.pRestart
                case 0
                    restartv = zerov;
                case 1
                    restartv = onev;
                otherwise
                    restartv = zerov;    % using default setting
            end
            for i = 2:(cCount.rw+1)
                lra(i,:) = lra(i-1,:)+ lr(i-1,:);      % adjusted or  restarted
                t(i,:)   = t(i-1,:) + onev;
                lcc = zerov;  % lower crossing count
                ucc = zerov;  % upper crossing count
                % lower crossing
                is_down = lra(i,:) < L;
                if sum(is_down) ~= 0
                    lcc(is_down) = floor(abs(lra(i,is_down)./ L(is_down)));      % lower crossing count
                end
                % upper crossing
                is_up = lra(i,:) > U;
                if sum(is_up) ~= 0
                    ucc(is_up)   = floor(abs(lra(i,is_up)  ./ U(is_up  )));      % upper crossing count
                end
                isbc = is_down + is_up;  % is boundary crossing
                if sum(isbc) ~= 0
                    % save timing of the events
                    isa = (isbc~=0);                   % select those columns where adustment is needed
                    cCount.pFET(i-1,isa) = t(i,isa);   % saving first exit time
                    t(i,isa)     = zerov(isa);         % setting back time to zero
                    % restarting process either at zero or at remainder
                    B            = lcc .* L + ucc .* U;    % value of boundary if relevant
                    lra(i,isa)   = (lra(i,isa) - B(isa)).* restartv(isa);     
                end
                % registering boundary crossing and convergence events
                cCount.pLcc(i, :) = lcc;
                cCount.pUcc(i, :) = ucc;
                cCount.pBCA(i, :) = cCount.pBCA(i-1, :) + ucc + lcc;
                cCount.pBCD(i, :) = cCount.pBCD(i-1, :) + ucc - lcc;
            end
            cCount.pTuc = sum(sum(cCount.pUcc));
            cCount.pTlc = sum(sum(cCount.pLcc));
            cCount.pTbc = cCount.pTuc +  cCount.pTlc;
            cCount.pP   = cCount.pTuc / cCount.pTbc;
            cCount.pQ   = 1 - cCount.pP; 
        end
        
        function mClassifyEvents(cCount)
            cCount.cTable      = zeros(2,3);
            cCount.cTable(1,1) = sum(sum(cCount.pLcc .* (cCount.pX <  0 ))) + 0.25;
            cCount.cTable(1,2) = sum(sum(cCount.pLcc .* (cCount.pX == 0 )));
            cCount.cTable(1,3) = sum(sum(cCount.pLcc .* (cCount.pX >  0 ))) + 0.25;
            cCount.cTable(2,1) = sum(sum(cCount.pUcc .* (cCount.pX <  0 ))) + 0.25;
            cCount.cTable(2,2) = sum(sum(cCount.pUcc .* (cCount.pX == 0 )));
            cCount.cTable(2,3) = sum(sum(cCount.pUcc .* (cCount.pX >  0 ))) + 0.25;
            outcomeSuccess     = cCount.cTable(1,3) + cCount.cTable(2,1); 
            outcomeTotal       = outcomeSuccess + cCount.cTable(1,1) + cCount.cTable(2,3);
            cCount.pPc         = outcomeSuccess / (outcomeSuccess +  outcomeTotal);
            cCount.pPVal       = p.s.binomTest(outcomeSuccess, outcomeTotal, 0.5 ,'Greater');
            cCount.pD          = cCount.pPVal < cCount.pAlpha;
        end
            
        %% unit root tests
        
        function mUnitRootFirstGenOptimalBoundary(cCount)
            cCount.mOptimalBoundary();
            cCount.mUnitRootFirstGen();
        end
        
        function mUnitRootFirstGen(cCount)
            cCount.mCount();
            cCount.pX           = zeros(cCount.rw+1, cCount.cl); 
            cCount.pX(2:end, :) = cCount.pBCD(1:end-1, :);
            cCount.mClassifyEvents(); 
        end
        
        % function mUnitRootSecondGen(cCount)
            % cCount.mOptimalBoundary();
            % cCount.mCount();
            % cCount.mClassifyEvents();
            %% toDo
            %% 'Not yet implemented'
        % end
        
        %% iterate functions
        
        function resTmp = mIterateOverBoundaries(cCount, cDgp, resTmp, i)
            if numel(cCount.pL.pOptions) ~= numel(cCount.pU.pOptions)
                error('Number of options for pL and pU must match');
            end
            cCount.pDX = cDgp.pDX;
            for j = 1:numel(cCount.pL.pOptions)
                cCount.pL.pV = cCount.pL.pOptions(j);
                cCount.pU.pV = cCount.pU.pOptions(j);
                cCount.mUnitRootFirstGen;
                resTmp(i, j) = cCount.pD;
            end
        end
            
    end
    
    methods(Access = private)
        
        function mSetDx(cCount)
            cCount.pStd = nanstd(cCount.pDX);
            [cCount.rw, cCount.cl] = size(cCount.pDX);
        end
        
    end
    
    methods(Static)
       
        function cCount = mCreateObj()
            cCount = p.bcc.CCount(); 
        end
        
        function cCount = mTestUnitRoot()
            cArx = p.dgp.CArx.mCreateObj();
            cArx.pAlpha.pV = 1;
            cArx.mSimulate();
            cCount = p.bcc.CCount(cArx.pDX);
            cCount.pL.pV = -3;
            cCount.pU.pV = - cCount.pL.pV;
            cCount.mUnitRootFirstGenOptimalBoundary();
        end
        
        function cCount = mTestStationary()
            cArx = p.dgp.CArx.mCreateObj();
            cArx.pAlpha.pV = 0.8;
            cArx.mSimulate();
            cCount = p.bcc.CCount(cArx.pDX);
            cCount.pL.pV = -3;
            cCount.pU.pV = - cCount.pL.pV;
            cCount.mUnitRootFirstGenOptimalBoundary();
        end
        
        function [cBcc, cDgp] = mTestIterateOverBoundaries()
            % create and set counting parameters
            cBcc                 = p.bcc.CCount;
            cBcc.pL.pOptions     = (-1:-1:-6);
            cBcc.pU.pOptions     = - cBcc.pL.pOptions;
            cBcc.pRestart        = 0;
            % create and set dgp parameters
            cDgp                 = p.dgp.CArx;
            cDgp.pStd            = 1;
            cDgp.pOverSample     = 200;
            cDgp.pAlpha.pOptions = [1, 0.9];
            cDgp.pN.pOptions     = 1;
            cDgp.pT.pOptions     = [100, 500, 1000];
            cDgp.pHeaderOutput   = {'b1', 'b2', 'b3', 'b4', 'b5', 'b6'};
            cDgp.mIterateNoDist(cBcc, 'mIterateOverBoundaries');
        end
        
    end
    
end