%% simulated data generating process
classdef CDgp < handle
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
     
    % data size
    properties
        pT          % # of time dimension 
        pN          % # of cross sections
        pSimCount   % number of simulations
        pOverSample % discard observations 1:pOverSample
    end
    
    % data
    properties
        pDistType
        pX
        pDX
        pY
        pDY
    end
    
    % IO
    properties
        pExportDir
        pFileName
        pFullFileName
    end
    
    % orchestration
    properties
        pHeaderOutput % headerInput consists of valOption type properties with options filled out
        pHowToAggregate
        pColumnWidth
    end
    
    properties(SetAccess = protected)
        pHeader   % pHeader = pHeaderInput + pHeaderOutput
        pRes      % cell, pRowCount x pInputCount + pOutoutCount
        pResT     % table
        pRowCount 
        pOutputCount
        rw
        cl
        pParamArray
    end
    
    methods(Abstract)
        mSimulate(cDgp)
    end
        
    methods
    
        function cDgp = CDgp()
            cDgp.pT = d.CValOptions;
            cDgp.pN = d.CValOptions;
            cDgp.pDistType = d.CValOptions({'n', 't'}); 
            cDgp.rw = 0;
            cDgp.cl = 0;
            % default parameters
            cDgp.pDistType.pV    = 'n';
            cDgp.pHowToAggregate = 'mAggregateByAverage';
            cDgp.pColumnWidth    = 35;
            cDgp.pParamArray     = cell(0,0);
            cDgp.pSimCount       = 200;
        end

        function mGetLevels(cDgp)
            cDgp.pX = cumsum(cDgp.pDX, 1);
        end
        
        function mGetFirstDifference(cDgp)
            if sum(sum(isnan(cDgp.pX)))==0
                cDgp.pDX = cDgp.pX(2:end,:) - cDgp.pX(1:(end-1),:); 
            else
                cDgp.pDX = a.NaNdif(cDgp.pDX);
            end
        end
        
        %% aggregate functions
        
        function mAggregateByAverage(cDgp, resTmp)
            tmp = mean(resTmp);
            for i = 1:numel(tmp)
                cDgp.pRes{cDgp.rw, cDgp.cl+i} = tmp(i);
            end 
        end
        
        %% orchestration
        
        function set.pHeaderOutput(cDgp, v)
            if isvector(v) == 0
                error(['Header for output must be a vector while u provided ', u.matrixSizeStr(v)]);
            end
            res = cell(2, numel(v));
            resT = cell2table(res);
            resT.Properties.VariableNames = v; % throws error is variable names are not ok
            cDgp.pHeaderOutput = v;
        end
        
        function mInitIterations(cDgp)
            if (cDgp.rw * cDgp.cl) == 0
                [v, ~, strCell] = mCountPropertiesWithOptions(cDgp);
                disp(['Iterating over ', u.cell2Str(strCell)]);
                disp(['Total casecount ', num2str(v)]);
                cDgp.rw = cDgp.rw + 1;
                cDgp.cl = cDgp.cl + 1;
                cDgp.pOutputCount = numel(cDgp.pHeaderOutput);
                if cDgp.pOutputCount < 1
                    error('At least 1 output is required. Pls fill out obj.pHeaderOutput.')
                end
                [rwTotal, clInput, ~] = mCountPropertiesWithOptions(cDgp);
                cDgp.pRes = cell(rwTotal, (clInput + cDgp.pOutputCount)); 
                if u.cl(cDgp.pHeaderOutput) ~= 1
                    cDgp.pHeaderOutput = cDgp.pHeaderOutput';
                end
            end  
        end
        
        function mAddParameter(cDgp, pStr, i)
            if iscell(cDgp.(pStr).pOptions)
                cDgp.(pStr).pV  = cDgp.(pStr).pOptions{i};
            else
                cDgp.(pStr).pV  = cDgp.(pStr).pOptions(i);
            end
            cDgp.cl = cDgp.cl + 1;
            if ismember(pStr, cDgp.pParamArray) == 0
                cc = numel(cDgp.pParamArray) + 1;
                cDgp.pParamArray{cc,1} = pStr;
            end
        end
        
        function mDoSimulations(cDgp, obj, methodStr)
            cDgp.mDisplayStep();
            if cDgp.rw == 1
                cDgp.pHeader = [cDgp.pParamArray; cDgp.pHeaderOutput];
            end
            % write input to results
            for i = 1:numel(cDgp.pParamArray)
                pStr = cDgp.pParamArray{i};
                if numel(cDgp.(pStr).pV) == 1
                    cDgp.pRes{cDgp.rw, i} = cDgp.(pStr).pV;
                else
                    if iscell(cDgp.(pStr).pV)
                        cDgp.pRes{cDgp.rw, i} = u.cell2Str(cDgp.(pStr).pV);
                    else
                        cDgp.pRes{cDgp.rw, i} = mat2str(cDgp.(pStr).pV);
                    end
                end
            end
            cDgp.cl = numel(cDgp.pParamArray);
            % simulate and write output to results
            resTmp = zeros(cDgp.pSimCount, cDgp.pOutputCount);
            for i = 1:cDgp.pSimCount
                cDgp.mSimulate();
                resTmp = obj.(methodStr)(cDgp, resTmp, i);
            end
            cDgp.(cDgp.pHowToAggregate)(resTmp);
            % set counters
            if cDgp.rw < cDgp.pRowCount
                cDgp.rw = cDgp.rw + 1;
                cDgp.cl = 1;
            else
                % setting back to the initial state
                cDgp.rw = 0;
                cDgp.cl = 0;
                cDgp.pParamArray = cell(0,0);
                cDgp.pResT = cell2table(cDgp.pRes);
                cDgp.pResT.Properties.VariableNames = cDgp.pHeader;
                disp(' ')
            end
        end
        
        function [v, cl, str] = mCountPropertiesWithOptions(cDgp)
            v = 1;
            cl = 0;
            str = cell(0, 0);
            tmp = properties(cDgp);
            for i = 1:numel(tmp)
                if ismember(tmp{i}, {'pRowCount', 'pClCount'}) == 0
                    if isa(cDgp.(tmp{i}), 'd.CValOptions')
                        if isempty(cDgp.(tmp{i}).pOptions) == 0
                            cl         = cl + 1;
                            v          = v * numel(cDgp.(tmp{i}).pOptions);
                            str{cl, 1} = tmp{i};
                        end
                    end
                end
            end
            if cl == 0
                error('None of the options are set. Pls. set at least one obj.property.options')
            end
            cDgp.pRowCount = v;
        end
        
        %% file IO
        
        function set.pExportDir(cDgp, v)
            if exist(v, 'file') ~= 7
                try mkdir(v)
                catch
                    error('Dir not exists and could not be created')
                end
            end
            cDgp.pExportDir = v;
        end
        
        function set.pFileName(cDgp, v)
            cDgp.pFileName = v;
            cDgp.mSetFullName();
        end
        
        function mSetFullName(cDgp)
            cDgp.pFullFileName = [cDgp.pExportDir, filesep, cDgp.pFileName];
        end
        
        function mExportRes(cDgp)
            writetable(cDgp.pResT, cDgp.pFullFileName);
        end 
            
    end
    
    methods(Access = private)
       
        % display result
        
        function mDisplayStep(cDgp)
            fprintf([num2str(cDgp.rw), ' ']);
            if cDgp.rw == cDgp.pColumnWidth
                disp(' ')
            end
        end
        
    end
    
end