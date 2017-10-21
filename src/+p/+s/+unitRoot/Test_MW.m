% PURPOSE: Maddala_Wu(1999) Test of Unit Root "A Comparative Study of Unit Root Tests with Panel Data and a New simple test"
% Oxford Bulletinb of Economics and Statistics, Special Issue, 631-652
%
% -------------------------------------------------------
% Usage:   MW(Y,Model,Lag_Orders,pmax)
% Where:   Y = matrix (T,N) of observations 
%          The data matrix can be not balanced. Missing Values must be specified as NaN
%
% Options
%       - Model = 1 : no individual effect 
%                 2 : individual effects (Default)
%                 3 : individual effects and time trends 
%       - Lag_Orders = if Lag_Orders is not specified the program determines
%                      the optimal lag order in individual ADF for each country (pmax=12 or T/4)
%                      In other cases, Lag_Orders must a vector (N,1) or (1,N) with
%                      the optimal lags for all the individuals of the panel                            
%       - pmax = Maximum Lag Order for individual ADF regressions 
%                     
% -------------------------------------------------------
% RETURNS:
% Panel unit root tests
%          MW.PMW                     Fisher statistic based on Individual ADF statistics
%          MW.PMW_Critical            Critical Values of the Fisher statistic at 1%, 5% and 10% 
%          MW.PMW_pvalue              Pvalue Pooled test statistic (Maddala Wu 1999)
%          MW.ZMW                     Choi (2001) statistic based on Individual ADF statistics (for large N tends to N(0,1)) 
%          MW.ZMW_Critical            Critical Values of the pooled test statistic (Choi, 2001) at 1%, 5% and 10% 
%          MW.ZMW_pvalue              Pvalue Pooled test statistic (Choi 2001)
%
% Individual ADF Tests
%          MW.pvalue                  Individual pvalues of individual ADF statitics
%          MW.pi                      Individual lag order in individual ADF models
%          MW.pmax                    Maximum Lag Order for individual ADF regressions 
%          MW.Ti                      Adjusted Individual Size
%          MW.tstats                  Individual ADF statistics
%          MW.sample                  Starting and Ending Dates of Adjusted Sample%           
%
% -------------------------------------------------------
%
% C. Hurlin, 01 June 2004
% LEO, University of Orléans
%

function [MW]=MW(Y,model,Lag_Orders,pmax);

if nargin==0; error (' Not enough argument for function MW(Y,model,Lag_Orders,pmax) ');end

if nargin==1; model=2; Lag_Orders=NaN; pmax=NaN ; end

if (model<1)|(model>3); model=2; end 

if nargin==2;  Lag_Orders=NaN ; pmax=NaN ; end

if nargin==3;  pmax=NaN ; end

if isnan(Lag_Orders)==0;  
    
    if length(Lag_Orders)~=size(Y,2) 

        warning ( ' The lag order must be specified for all the individuals of the panel : your lag structure cannot be considered ')
        
        Lag_Orders=NaN;
        
    end
    
end

%------------------------
%--- Transformed Data ---
%------------------------

T=size(Y,1);                                               % Time Dimension 

N=size(Y,2);                                               % Individual Dimension

if isnan(Lag_Orders)
    
    ADF = p.s.unitRoot.ADF_Individual(Y,model,NaN,pmax);                  % Individual ADF REgression for get Lag Order ADF.pi and size ADF.Ti
    
else 
    ADF = p.s.unitRoot.ADF_Individual(Y,model,Lag_Orders);                % Individual ADF REgression (the lag structure is specified)
    
end
    
dy=Y(2:end,:)-Y(1:end-1,:);                                % Vector of First Order Differences     

DYlag=ones(T,N,max(ADF.pi))*NaN;                           % Matrix of Lagged First Order Differences     

for p=1:max(ADF.pi)           

    DYlag(p+2:end,:,p)=dy(1:end-p,:);                      % Matrix of Lagged First Order Differences       
    
end

%----------------------------------------------------------
%--- Test 1 : Fisher Based on ADF individual Statistics ---
%----------------------------------------------------------
MW.PMW=-2*sum(log(ADF.pvalue));                            % Fisher statistic based on Individual ADF statistics

MW.PMW_Critical=chi2inv([0.99 0.95 0.90],2*N);            % Critical Values of the Fisher statistic at 1%, 5% and 10% 

MW.PMW_pvalue= 1-gammainc(MW.PMW/2, N);                    % Pvalue Pooled test statistic (Maddala Wu 1999)

MW.ZMW=(1/(2*sqrt(N)))*sum(-2*log(ADF.pvalue)-2);          % Choi (2001) statistic based on Individual ADF statistics (for large N tends to N(0,1)) 

MW.ZMW_Critical=[2.3263 1.6449 1.2816];                    % Critical Values of the pooled test statistic (Choi, 2001) at 1%, 5% and 10% 

MW.ZMW_pvalue=1-normcdf(MW.ZMW);                          % Pvalue Pooled test statistic (Choi 2001)

%===============
%===============
%=== RESULTS ===
%===============
%===============

MW.pvalue=ADF.pvalue;                                      % Individual pvalues of individual ADF statitics

if isnan(Lag_Orders)
    
    MW.pi=ADF.pi;                                          % Individual lag order in individual ADF models
    
else

    MW.pi=Lag_Orders(:);                                   % Individual lag order in individual ADF models    
    
end  

MW.pmax=ADF.pmax;                                          % Maximum of individual lag order in individual ADF models    

MW.Ti=ADF.Ti;                                              % Adjusted Individual Size

MW.tstats=ADF.tstat;                                       % Individual ADF statistics

MW.sample=[ADF.ai ADF.bi];                                 % Starting and Ending Dates of Adjusted Sample

switch model
    
    case 1; MW.model=' Model without intercept nor trend: Model 1';
    
    case 2; MW.model=' Model with intercept: Model 2';
    
    case 3; MW.model=' Model with intercept and trend: Model 3 ';
    
end

%
% End of Program
%

