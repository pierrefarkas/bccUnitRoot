% PURPOSE: Im, Pesaran and Shin (2003)test of unit root
% "Testing for Unit Root in Heterogeneous Panels", Journal of Econometrics, 115, 53-74
% 3.
% -------------------------------------------------------
% Usage:   Test_IPS(Y,Model,pmax)
% Where:   Y = matrix (T,N) of observations 
%          The data matrix can be not balanced. Missing Values must be specified as NaN 
% Options
%       - Model = 1 : no individual effect 
%                 2 : individual effects (Default)
%                 3 : individual effects and time trends 
%       - pmax = Maximum Lag Order for individual ADF regressions 
%
% -------------------------------------------------------
% RETURNS:
%          results.tbar                    Mean of Individual Augmented Dickey Fuller statistics              
%          results.Wbar                    Standardized IPS statistic based on E[ti(pi,0)/bi=0] and V[ti(pi,0)/bi=0]
%          results.Wbar_pvalue             Pvalue of Wbar
%          results.Zbar                    Standardized IPS statistic based on the moments of the DF distribution 
%          results.Zbar_pvalue             Pvalue of Zbar
%          results.critical                Critical Values of the Normal distribution at 1%, 5% and 10%
%          results.tbar_DF                 Mean of Individual Dickey Fuller statistics              
%          results.Zbar_DF                 Standardized IPS statistic (assumption of no autocorrelation of residuals)
%          results.Zbar_DF_pvalue          Pvalue of Zbar_DF
%          results.pi                      Individual lag order in individual ADF models
%          results.pmax                    Maximum of lag order in individual ADF models
%          results.Ti                      Adjusted Individual Size
%           
% -------------------------------------------------------
%
% C. Hurlin, 25 Mai 2004
% LEO, University of Orléans
%

function [results]=Test_IPS(Y,model,pmax);

if nargin==0; error (' Not enough argument for function Test_IPS(Y,model,pmax) ');end

if nargin==1; model=2; pmax=NaN ;end

if nargin==2; model=2; pmax=NaN ;end

if (model<1)|(model>3);model=2;end 

%------------------------
%--- Transformed Data ---
%------------------------

T=size(Y,1);                           % Time Dimension 

N=size(Y,2);                           % Individual Dimension

ADF = p.s.unitRoot.ADF_Individual(Y,model,NaN,pmax);  % Individual ADF REgression for get Lag Order ADF.pi and size ADF.Ti

dy=Y(2:end,:)-Y(1:end-1,:);            % Vector of First Order Differences     

DYlag=ones(T,N,max(ADF.pi))*NaN;       % Matrix of Lagged First Order Differences     

for p=1:max(ADF.pi)           

    DYlag(p+2:end,:,p)=dy(1:end-p,:);  % Matrix of Lagged First Order Differences       
    
end

dy=[NaN*ones(1,N) ; dy];                % Vector of First Order Differences     


%==============================================================================
%==============================================================================
%=== PART I : Test under the assumption of no Auto-Correlation of Residuals ===
%==============================================================================
%==============================================================================

%--------------------------------
%--- Individual DF Statistics ---
%--------------------------------

DF=zeros(N,1);                                                                        % Individual DF Statistics

for i=1:N

    sample=[dy(:,i) [NaN ; Y(1:end-1,i)] ];                                           % Baseline Data (Selection of Sample)

    SAMPLE=zeros(size(sample));SAMPLE(isnan(sample)==1)=1;
    
    ai=max(sum(cumprod(SAMPLE)))+1;                                                   % Starting Date of Adjusted Sample

    [a,b]=max(cumsum(ones(size(sample))-isnan(sample)));bi=min(b);                    % Ending Date of Adjusted Sa
            
    Tip=bi-ai+1;                                                                      % Adjusted Size for Individual i

    if Tip<10 
        
        disp(sprintf('  WARNING : Individual time dimensions must be superior to 10 (or 6 for balanced data) : see unit n°%1g ',i))
        
    end 
        
    Xi=Y(ai-1:bi-1,i);                                                                % Matrix of Regressors
       
    switch model              
        
        case 2; Xi=[ones(Tip,1) Xi];                                                  % 'Model with intercept'
            
        case 3; Xi=[ones(Tip,1) (1:1:Tip)' Xi];                                       % 'Model with intercept and trend'
     
    end
    
    if sum(isnan( Xi(:,2) ) + isnan(dy(ai:bi,i))) > 0
        isa = isnan( Xi(:,2) ) + isnan(dy(ai:bi,i))==0;
        dy0 = dy(ai:bi,i);
        dy1 = dy0(isa);
        Xi1 = Xi(isa,:);
        coef_ind    = Xi1\dy1;                                                           % Coefficients ADF Regression     
        var_res_ind = sum((dy1-Xi1*coef_ind).^2)/(Tip-model);                         % Residual Variance of ADF
        tstat_ind   = coef_ind./diag(sqrt(var_res_ind*inv(Xi1'*Xi1)));                  % Tstat ADF Regression
    else
       % no missing values
        coef_ind    = Xi\dy(ai:bi,i);                                                       % Coefficients DF Regression     
        var_res_ind = sum((dy(ai:bi,i)-Xi*coef_ind).^2)/(Tip-model);                     % Residual Variance of DF
        tstat_ind   = coef_ind./diag(sqrt(var_res_ind*inv(Xi'*Xi)));                       % Tstat DF Regression            
    end
    DF(i,1)=tstat_ind(end,1);                                                          % Individual DF statistic
    T_DF(i,1)=Tip;                                                                     % Adjusted sample size for DF model
    
end        
   
%--------------------------------------------------------------------------------
%--- Moments of Individual ADF statistics under assumption of no correlation ----
%--------------------------------------------------------------------------------
% The statistic used here is the standard t(it) and not t(it)tild (page 4). 
% The corresponding moments from table 1, page 60, are used. 

TT_DF=[6 7 8 9 10 15 20 25 30 40 50 100 500 1000 5000]';

esperance=[-1.520 -1.515 -1.501 -1.501 -1.504 -1.514 -1.522 -1.520 -1.526 -1.523 -1.527 -1.532 -1.531 -1.529 -1.533]';

variance=[1.745 1.414 1.228 1.132 1.069 0.923 0.851 0.809 0.789 0.770 0.760 0.735 0.715 0.707 0.706]';

Em_DF=0;                                                            % Mean of Individual Esperances of DF tstats 

Vm_DF=0;                                                            % Mean of Individual Variances of DF tstats

for i=1:N                                                           % Computation of the Moments of the individual DF statistic

    [a,posi_DF]=min((TT_DF-T_DF(i,1)).^2);                          % Posi : Position in TT table
    
    Em_DF=Em_DF+esperance(posi_DF);                                 % Sum of Individual Esperances of DF tstats 
    
    Vm_DF=Vm_DF+variance(posi_DF);                                  % Sum of Individual Variances of DF tstats 
    
end

%---------------------------
%--- Statistics of Tests ---
%---------------------------
tbar_DF=mean(DF);                                                    % Mean of Individual Dickey Fuller statistics

Zbar_DF=sqrt(N)*(tbar_DF-Em_DF/N)/sqrt(Vm_DF/N);                     % Standardized IPS statistic (case I : no autocorrelation)


%=============================================================================
%=============================================================================
%=== PART II : Test in the general case with Auto-Correlation of Residuals ===
%=============================================================================
%=============================================================================

%---------------------------------------------
%--- Moments of Individual ADF statistics ----
%---------------------------------------------
% The corresponding moments are issued from table 2, page 18

TT=[10 15 20 25 30 40 50 60 70 100]';

PP=(0:1:8)';

m2=[-1.504 -1.514 -1.522 -1.520 -1.526 -1.523 -1.527 -1.519 -1.524 -1.532 ;
    -1.488 -1.503 -1.516 -1.514 -1.519 -1.520 -1.524 -1.519 -1.522 -1.530 ;
    -1.319 -1.387 -1.428 -1.443 -1.460 -1.476 -1.493 -1.490 -1.498 -1.514 ;
    -1.306 -1.366 -1.413 -1.433 -1.453 -1.471 -1.489 -1.486 -1.495 -1.512 ;
    -1.171 -1.260 -1.329 -1.363 -1.394 -1.428 -1.454 -1.458 -1.470 -1.495 ;
       NaN    NaN -1.313 -1.351 -1.384 -1.421 -1.451 -1.454 -1.467 -1.494 ;
       NaN    NaN    NaN -1.289 -1.331 -1.380 -1.418 -1.427 -1.444 -1.476 ;
       NaN    NaN    NaN -1.273 -1.319 -1.371 -1.411 -1.423 -1.441 -1.474 ;
       NaN    NaN    NaN -1.212 -1.266 -1.329 -1.377 -1.393 -1.415 -1.456 ] ;
       
v2=[1.069 0.923 0.851 0.809 0.789 0.770 0.760 0.749 0.736 0.735 ;
    1.255 1.011 0.915 0.861 0.831 0.803 0.781 0.770 0.753 0.745 ;
    1.421 1.078 0.969 0.905 0.865 0.830 0.798 0.789 0.766 0.754 ;
    1.759 1.181 1.037 0.952 0.907 0.858 0.819 0.802 0.782 0.761 ;
    2.080 1.279 1.097 1.005 0.946 0.886 0.842 0.819 0.801 0.771 ;
      NaN   NaN 1.171 1.055 0.980 0.912 0.863 0.839 0.814 0.781 ;
      NaN   NaN   NaN 1.114 1.023 0.942 0.886 0.858 0.834 0.795 ;
      NaN   NaN   NaN 1.164 1.062 0.968 0.910 0.875 0.851 0.806 ;    
      NaN   NaN   NaN 1.217 1.105 0.996 0.929 0.896 0.871 0.818 ];
 
m3=[-2.166 -2.167 -2.168 -2.167 -2.172 -2.173 -2.176 -2.174 -2.174 -2.177 ;
    -2.173 -2.169 -2.172 -2.172 -2.173 -2.177 -2.180 -2.178 -2.176 -2.179 ;
    -1.914 -1.999 -2.047 -2.074 -2.095 -2.120 -2.137 -2.143 -2.146 -2.158 ;
    -1.922 -1.977 -2.032 -2.065 -2.091 -2.117 -2.137 -2.142 -2.146 -2.158 ;
    -1.750 -1.823 -1.911 -1.968 -2.009 -2.057 -2.091 -2.103 -2.114 -2.135 ;
       NaN    NaN -1.888 -1.955 -1.998 -2.051 -2.087 -2.101 -2.111 -2.135 ;
       NaN    NaN    NaN -1.868 -1.923 -1.995 -2.042 -2.065 -2.081 -2.113 ;
       NaN    NaN    NaN -1.851 -1.912 -1.986 -2.036 -2.063 -2.079 -2.112 ;
       NaN    NaN    NaN -1.761 -1.835 -1.925 -1.987 -2.024 -2.046 -2.088 ];
      
v3=[1.132 0.869 0.763 0.713 0.690 0.655 0.633 0.621 0.610 0.597 ;
    1.453 0.975 0.845 0.769 0.734 0.687 0.654 0.641 0.627 0.605 ;
    1.627 1.036 0.882 0.796 0.756 0.702 0.661 0.653 0.634 0.613 ;
    2.482 1.214 0.983 0.861 0.808 0.735 0.688 0.674 0.650 0.625 ;
    3.947 1.332 1.052 0.913 0.845 0.759 0.705 0.685 0.662 0.629 ;
      NaN   NaN 1.165 0.991 0.899 0.792 0.730 0.705 0.673 0.638 ;
      NaN   NaN   NaN 1.055 0.945 0.828 0.753 0.725 0.689 0.650 ;
      NaN   NaN   NaN 1.145 1.009 0.872 0.786 0.747 0.713 0.661 ;
      NaN   NaN   NaN 1.208 1.063 0.902 0.808 0.766 0.728 0.670];

%-------------------------------------------------------------
%--- Standardisation E[ti(pi,0)/bi=0] and V[ti(pi,0)/bi=0] ---
%-------------------------------------------------------------

Em=0;                                                               % Mean of Individual Esperances E[ti(pi,0)/bi=0]

Vm=0;                                                               % Mean of Individual Variances V[ti(pi,0)/bi=0]

for i=1:N                                                           % Computation of the Moments of the individual ADF statistic

    [a,T_posi]=min((TT-ADF.Ti(i,1)).^2);                            % Posi : Position in TT table
    
    [b,p_posi]=min((PP-ADF.pi(i,1)).^2);                            % Posi : Position in TT table

    switch model
        
    case 1, Em=NaN; Vm=NaN;                                          % Case not tabulated in (IPS 2003)        
         
        disp('  WARNING : the moments E[ti(pi,0)/bi=0] and V[ti(pi,0)/bi=0] are not tabulated for Model 1 in IPS')
        
    case 2

        Em=Em+m2(p_posi,T_posi);                                     % Sum of Individual Esperances of DF tstats             

        Vm=Vm+v2(p_posi,T_posi);                                     % Sum of Individual Esperances of DF tstats             
       
        if isnan(m2(p_posi,T_posi)+v2(p_posi,T_posi))==1
        
            disp(sprintf('  WARNING : Given the optimal lag structure, individual time dimension is not sufficient for unit n°%1g : T = %1g for p = %1g ',i,ADF.Ti(i),ADF.pi(i)))
        
        end 
        
    case 3

        Em=Em+m3(p_posi,T_posi);                                     % Sum of Individual Esperances of DF tstats             

        Vm=Vm+v3(p_posi,T_posi);                                     % Sum of Individual Esperances of DF tstats             

        if isnan(m3(p_posi,T_posi)+v3(p_posi,T_posi))==1
        
            disp(sprintf('  WARNING : Given the optimal lag structure, individual time dimension is not sufficient for unit n°%1g : T = %1g for p = %1g ',i,ADF.Ti(i),ADF.pi(i)))
        
        end 
        
    end
    
end
      
      
%------------------------------------------
%--- Statistics of Test : tbar and Wbar ---
%------------------------------------------

tbar=mean(ADF.tstat);                                              % Mean of Individual Augmented Dickey Fuller statistics
      
Zbar=sqrt(N)*(tbar-Em_DF/N)/sqrt(Vm_DF/N);                         % Standardized IPS statistic 

Wbar=sqrt(N)*(tbar-Em/N)/sqrt(Vm/N);                               % Standardized IPS statistic based on E[ti(pi,0)/bi=0] and V[ti(pi,0)/bi=0]


%===============
%===============
%=== RESULTS ===
%===============
%===============

results.tbar=tbar;                                              % Mean of Individual Augmented Dickey Fuller statistics              

results.Wbar=Wbar;                                              % Standardized IPS statistic based on E[ti(pi,0)/bi=0] and V[ti(pi,0)/bi=0]

asdfasdf = 0;

results.Wbar_pvalue=normcdf(Wbar);                             % Pvalue of Wbar

results.Zbar=Zbar;                                              % Standardized IPS statistic based on the moments of the DF distribution 

results.Zbar_pvalue=normcdf(Zbar);                             % Pvalue of Zbar

results.critical= -[2.3263 1.6449 1.2816 ];                     % Critical Values of the Normal distribution at 1%, 5% and 10%

results.tbar_DF=tbar_DF;                                        % Mean of Individual Dickey Fuller statistics              

results.Zbar_DF=Zbar_DF;                                        % Standardized IPS statistic (assumption of no autocorrelation of residuals)

results.Zbar_DF_pvalue=normcdf(Zbar_DF);                       % Pvalue of Zbar_DF

results.pi=ADF.pi;                                              % Individual lag order in individual ADF models

results.tstats_ind=ADF.tstat;                                   % Individual lag order in individual ADF models

results.pmax=pmax;                                              % Maximum of lag order in individual ADF models

results.Ti=ADF.Ti;                                              % Adjusted Individual Size

switch model
    
    case 1; results.model=' Model without intercept nor trend: Model 1';
    
    case 2; results.model=' Model with intercept: Model 2';
    
    case 3; results.model=' Model with intercept and trend: Model 3 ';
    
end

%
% End of Program
%

