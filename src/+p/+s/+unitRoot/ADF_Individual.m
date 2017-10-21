% PURPOSE: ADF (1981) Unit Root Tests on Individual Times Series
% Augmented Dickey Fuller Tests : Lag Selection with the k-max Criterium
% Dickey et Fuller (1981), Econometrica
%
% -------------------------------------------------------
% Usage:   ADF_Individual(Y,Model,Lag_order,pmax)
%
% Where:   Y = matrix (T,N) of observations 
%          The data matrix can be not balanced. Missing Values must be specified as NaN
%
% Options
%       - Model = 1 : no individual effect 
%                 2 : individual effects (Default)
%                 3 : individual effects and time trends 
%       - Lag_Orders = if Lag_Orders is not specified the program determines
%                      the optimal lag order for each country (pmax=12 or T/4)
%                      In other cases, Lag_Orders must a vector (N,1) or (1,N) with
%                      the optimal lags for all the individuals of the panel
%                      The criteria used here is the Schwartz criteria
%       - pmax = Maximum of the lag order authorized
%                              
% -------------------------------------------------------
% RETURNS:
%           ADF.tstat           Individal ADF Statistic
%           ADF.pvalue          Pvalue for each individual ADF Statistic
%           ADF.critical        Critical Values of the DF distribution at 1%, 5% and 10% for T and N sample
%           ADF.pi              Optimal Lag Order in ADF Regression
%           ADF.Ti              Optimal Adjusted Size (Unbalanced Panel)
%           ADF.rho             Estimated autoregressive parameter
%           ADF.sig             Standard Error of ADF Residual
%           ADF.residual        ADF Residual
%           ADF.model           Model of the ADF test
%           ADF.si;             Individual Ratios of long run standard deviation to the innovation SE
%           ADF.ai              Starting Date of Adjusted Sample
%           ADF.bi              Ending Date of Adjusted Sample    
%           ADF.pmax            Maximum of the lag order authorized
%
% Remark : This program can consider different samples for each individual (non balanced panel)
%
% -------------------------------------------------------
%
% C. Hurlin, June 2006
% LEO, University of Orléans
%

function [ADF]=ADF_Individual(Y,model,Lag_Orders,pmax);

% disp('do I modify the correct thing?')
% asdfasdfa = 0;    

% persistent dfstat_p;
strPath = mfilename('fullpath');
dfstat_p = load([strPath, filesep, 'Critical_DF.mat']);

if nargin==1; model=2; Lag_Orders=NaN; pmax=NaN; end

if nargin==2;  Lag_Orders=NaN ; pmax=NaN; end

if nargin==3;  pmax=NaN; end

if (model<1)|(model>3);model=2;end 

if isnan(Lag_Orders)==0;  


    
    if length(Lag_Orders)~=size(Y,2) 

        warning ( ' The lag order must be specified for all the individuals of the panel : your lag structure cannot be considered ')
        
        Lag_Orders=NaN;
        
    end
    
end


%------------------------
%--- Transformed Data ---
%------------------------

T=size(Y,1);                                                                        % Time Dimension 

N=size(Y,2);                                                                        % Individual Dimension

if isnan(pmax)==1
    
    pmax=min(floor(T/4),12);                                                        % Maximum of Lag order
    
end

dy=Y(2:end,:)-Y(1:end-1,:);                                                         % Vector of First Order Differences     

DYlag=ones(T,N,pmax)*NaN;                                                           % Matrix of Lagged First Order Differences     

for p=1:pmax                                                                        % Loop on lag order

    DYlag(p+2:end,:,p)=dy(1:end-p,:);                                               % Matrix of Lagged First Order Differences     
    
end                                                                                 % End of loop


%------------------------
%--- Data for Pvalues ---
%------------------------

switch model                                                                        % Switch on deterministic component 
    
    case 1, critical = dfstat_p.Critical_M1;                                                   % Critical values model 1
            
    case 2, critical = dfstat_p.Critical_M2;                                                   % Critical values model 2
        
    case 3, critical = dfstat_p.Critical_M3;                                                   % Critical values model 3
        
end                                                                                 % End of switch 

TT_pvalue=critical(1,2:end)';                                                       % Sample dimension in tables

proba=critical(2:end,1)';                                                           % Probabilities used in tables

critical=critical(2:end,2:end);                                                     % Critical Values

%---------------------------------------
%--- Individual Lag Length Selection ---
%---------------------------------------

ADF.tstat=zeros(N,1);                                                               % Vector of Estimated ADF stat

ADF.critical=zeros(N,3);                                                            % Critical Values of the DF distribution at 1%, 5% and 10% for T and N sample

ADF.pvalue=zeros(N,1);                                                              % Vector of Pvalues

ADF.pi=zeros(N,1);                                                                  % Vector of Individual Lag Length

ADF.Ti=zeros(N,1);                                                                  % Vector of Individual Adjusted Size

ADF.rho=zeros(N,1);                                                                 % Vector of estimated autoregressive parameter

ADF.sig=zeros(N,1);                                                                 % Vector of Standard Error of ADF Residual

ADF.residual=ones(T,N)*NaN;                                                         % ADF Residual

ADF.si=zeros(N,1);                                                                  % Individual Ratios of long run standard deviation to the innovation SE

ADF.ai=zeros(N,1);                                                                  % Starting Date of Adjusted Sample

ADF.bi=zeros(N,1);                                                                  % Ending Date of Adjusted Sample    

dy=[NaN*ones(1,N) ; dy];                                                            % Vector of First Order Differences     

for i=1:N                                                                           % Loop on individual unit

    p_indi=0;                                                                       % Initial Guess on the Optimal Lag order 
   
    tstat_ind=zeros(p_indi+2,1);                                                    % Initialisation of the vector of tstats                                  
    
    base=[dy(:,i) [NaN ; Y(1:end-1,i)] ];                                           % Baseline Data (Selection of Sample)
        
    Tip=T-1;                                                                        % Adjusted T size (initialisation)

    if prod(double(isnan(Lag_Orders)==0))==1

        finboucle=1;                                                                % Particular case : the user gives the lag structure
        
    else
            
        finboucle=pmax+1;                                                           % Gerneral case : the user do not give the lag structure
    
    end

    for indic_j=1:finboucle                                                         % Criteria Scwartz Information Criteria

        p_indi=indic_j-1;                                                           % Lag Order
        
        if isnan(Lag_Orders)==0 ; p_indi=Lag_Orders(i); ; end                       % Particular case : the user gives the lag structure
                    
        if p_indi>0                                                                 % ADF case
            
            sample=[base reshape(DYlag(:,i,1:p_indi),T,p_indi)];                    % Initial Sample (logical variable)                      

            SAMPLE=zeros(size(sample));SAMPLE(isnan(sample)==1)=1;
            
            ai=max(sum(cumprod(SAMPLE)))+1;                                         % Starting Date of Adjusted Sample

            [a,b]=max(cumsum(ones(size(sample))-isnan(sample)));bi=min(b);          % Ending Date of Adjusted Sample

            Tip=bi-ai+1;                                                            % Adjusted Size for Individual i
                       
            Xi=[Y(ai-1:bi-1,i) reshape(DYlag(ai:bi,i,1:p_indi),Tip,p_indi)];        % Matrix of Regressors

        else                                                                        % DF case

            sample=base;                                                            % Initial Sample            
            
            SAMPLE=zeros(size(sample));SAMPLE(isnan(sample)==1)=1;
                      
            ai=max(sum(cumprod(SAMPLE)))+1;                                         % Starting Date of Adjusted Sample

            [a,b]=max(cumsum(ones(size(sample))-isnan(sample)));bi=min(b);          % Ending Date of Adjusted Sa
            
            Tip=bi-ai+1;                                                            % Adjusted Size for Individual i
        
            Xi=Y(ai-1:bi-1,i);                                                      % Matrix of Regressors
            
        end                                                                         % End of the loop

        %-----------------------------------------------
        %--- Individual Regression with pi ADF terms ---
        %-----------------------------------------------
    
        switch model                                                                % Switch on deterministic component
            
        case 2                                                                      % Model with intercept
            % disp('here')
            Xi=[ones(Tip,1) Xi];                                                    % Matrix of regressors
            
        case 3                                                                      % Model with intercept and trend
            
            Xi=[ones(Tip,1) (1:1:Tip)' Xi];                                         % Matrix of regressors
            
        end                                                                         % End of switch
        
        % asfdasfsa = 0; % entry point for debug: missing data
        
        % remove columns with missing values
        if sum(isnan( Xi(:,2) ) + isnan(dy(ai:bi,i))) > 0
            isa = isnan( Xi(:,2) ) + isnan(dy(ai:bi,i))==0;
            dy0 = dy(ai:bi,i);
            dy1 = dy0(isa);
            Xi1 = Xi(isa,:);
            coef_ind = Xi1\dy1;                                                           % Coefficients ADF Regression     
            var_res_ind = sum((dy1-Xi1*coef_ind).^2)/(Tip - p_indi-model);                % Residual Variance of ADF
            tstat_ind = coef_ind./diag(sqrt(var_res_ind * inv(Xi1'*Xi1)));                  % Tstat ADF Regression
            rss = sum((dy1- Xi1*coef_ind).^2);                                            % RSS
        else
            coef_ind = Xi\dy(ai:bi,i);                                                    % Coefficients ADF Regression     
            var_res_ind = sum((dy(ai:bi,i)-Xi*coef_ind).^2)/(Tip-p_indi-model);           % Residual Variance of ADF
            % if ai == 12
            %    disp('debug')
            % end
            tstat_ind = coef_ind./diag(sqrt(var_res_ind*inv(Xi'*Xi)));                      % Tstat ADF Regression            
            rss=sum((dy(ai:bi,i)-Xi*coef_ind).^2);                                        % RSS
        end
          
        nbparam=p_indi+model-1;                                                     % Number of Parameters

        AIC(indic_j,i)=log(rss/(Tip-nbparam-1))+2*(nbparam)/(Tip);                  % AIC Criteria

        BIC(indic_j,i)=log(rss/(Tip-nbparam-1))+(nbparam)/(Tip)*log(Tip);           % BIC Criteria

    end                                                                             % End of while
    
    %---------------------------------------
    %--- ADF Regression with Optimal Lag ---
    %---------------------------------------
    [a,b]=min(BIC(:,i));p_indi=b-1;                                                 % Optimal Lag
    
    if isnan(Lag_Orders)==0 ; p_indi=Lag_Orders(i); ; end                           % Particular case : the user gives the lag structure

    if p_indi>0                                                                     % ADF case
            
        sample=[base reshape(DYlag(:,i,1:p_indi),T,p_indi)];                        % Initial Sample (logical variable)                      

        SAMPLE=zeros(size(sample));SAMPLE(isnan(sample)==1)=1;
            
        ai=max(sum(cumprod(SAMPLE)))+1;                                             % Starting Date of Adjusted Sample

        [a,b]=max(cumsum(ones(size(sample))-isnan(sample)));bi=min(b);              % Ending Date of Adjusted Sample

        Tip=bi-ai+1;                                                                % Adjusted Size for Individual i
                       
        Xi=[Y(ai-1:bi-1,i) reshape(DYlag(ai:bi,i,1:p_indi),Tip,p_indi)];            % Matrix of Regressors

    else                                                                            % DF case

        sample=base;                                                                % Initial Sample            
            
        SAMPLE=zeros(size(sample));SAMPLE(isnan(sample)==1)=1;
                      
        ai=max(sum(cumprod(SAMPLE)))+1;                                             % Starting Date of Adjusted Sample

        [a,b]=max(cumsum(ones(size(sample))-isnan(sample)));bi=min(b);              % Ending Date of Adjusted Sa
            
        Tip=bi-ai+1;                                                                % Adjusted Size for Individual i
        
        Xi=Y(ai-1:bi-1,i);                                                          % Matrix of Regressors
            
    end                                                                             % End of the loop

    switch model                                                                    % Switch on deterministic component
            
    case 2                                                                          % Model with intercept
        
        Xi=[ones(Tip,1) Xi];                                                        % Matrix of regressors
            
    case 3                                                                          % Model with intercept and trend
            
        Xi=[ones(Tip,1) (1:1:Tip)' Xi];                                             % Matrix of regressors
            
    end                                                                             % End of switch
        
    if sum(isnan( Xi(:,2) ) + isnan(dy(ai:bi,i))) > 0
        isa = isnan( Xi(:,2) ) + isnan(dy(ai:bi,i))==0;
        dy0 = dy(ai:bi,i);
        dy1 = dy0(isa);
        Xi1 = Xi(isa,:);
        coef_ind = Xi1\dy1;                                                           % Coefficients ADF Regression     
        var_res_ind = sum((dy1 - Xi1 * coef_ind).^2)/(Tip - p_indi-model);                % Residual Variance of ADF
        tstat_ind = coef_ind./diag(sqrt(var_res_ind*inv(Xi1'*Xi1)));                  % Tstat ADF Regression
    else
        coef_ind=Xi\dy(ai:bi,i);                                                        % Coefficients ADF Regression     
        var_res_ind=sum((dy(ai:bi,i)-Xi*coef_ind).^2)/(Tip-p_indi-model);               % Residual Variance of ADF
        tstat_ind=coef_ind./diag(sqrt(var_res_ind*inv(Xi'*Xi)));                        % Tstat ADF Regression            
    end
        
        
    %--------------------------------
    %--- Pvalue of ADF statistics ---
    %--------------------------------
    
    [a,posi_T]=min((Tip-TT_pvalue).^2);                                             % Approximated sample size
        
    [a,posi_p]=min((tstat_ind(model)-critical(:,posi_T)).^2);                       % Approximated critical value
    
    pvalue=proba(posi_p);                                                           % Pvalue for unit i
            
    [d,posi_prob1]=min((proba-0.01).^2);                                            % Position of the approximated critical value at 1%
 
    [d,posi_prob5]=min((proba-0.05).^2);                                            % Position of the approximated critical value at 5%

    [d,posi_prob10]=min((proba-0.10).^2);                                           % Position of the approximated critical value at 10%

    %-----------------------------------------
    %--- Optimal Lag Order and ADF Results ---
    %-----------------------------------------

    ADF.tstat(i)=tstat_ind(model);                                                  % Individual ADF Statistic

    ADF.pvalue(i)=pvalue;                                                           % Pvalue for each individual ADF Statistic    
    
    switch model                                                                    % Switch on deterministic component
        
    case 1                                                                          % Model with no intercept
    
        ADF.model='Model without intercept nor trend';                              % Model with no intercept
        
    case 2                                                                          % Model with intercept
     
        ADF.model='Model with intercept';                                           % Model with intercept 
        
    case 3                                                                          % Model with intercept and trend
        
        ADF.model='Model with intercept and trend';                                 % Model with intercept and trend
        
    end                                                                             % End of switch
    
    if isnan(Lag_Orders)                                                            % Case of optimal lag structure
        
        ADF.pi(i)=p_indi;                                                           % Optimal Lag Order in ADF Regression
        
    else                                                                            % The structure lag is given by user
        
        ADF.pi(i)=Lag_Orders(i);                                                    % Optimal Lag Order in ADF Regression
        
    end                                                                             % End of if condition

    ADF.rho(i)=coef_ind(model);                                                     % Estimated autoregressive parameter
    
    ADF.Ti(i)=Tip;                                                                  % Optimal Adjusted Size (Unbalanced Panel)
    
    ADF.sig(i)=sqrt(var_res_ind);                                                   % Standard Error of ADF Individual Residuals

    ADF.residual(end-Tip+1:end,i)=dy(ai:bi,i)-Xi*coef_ind;                          % ADF Residual             

    ADF.si(i)=abs(1-sum(coef_ind(model+1:end)));                                    % Individual Ratios of long run standard deviation to the innovation SE    

    ADF.ai(i)=ai;                                                                   % Starting Date of Adjusted Sample

    ADF.bi(i)=bi;                                                                   % Ending Date of Adjusted Sample    
    
    ADF.critical(i,:)=critical([posi_prob1 posi_prob5 posi_prob10],posi_T)';        % Critical Values of the DF distribution at 1%, 5% and 10% for T and N sample
    
    ADF.pmax=pmax;                                                                  % Maximum of the lag order authorized
    
end                                                                                 % End of loop on individual

%
% End of Program
%
