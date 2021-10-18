%Main impulse-matching file
%
%Jarocinski-Karadi, 2019

clear all; close all; clc;
warning('off','all');

global oo_ M_ options_;

MainDir = pwd;
addpath(MainDir);
addpath([MainDir '/auxiliaryCode']);

%Switches for impulse response runs. 
%[For optimization, 
%           First run with 'CHOL' and 'MonPol'
%           Then with 'BL' and 'MonPol' and flagKappa0 = 0
%           Then change starting value for flagKappa0=1 (line 172) and run with flagKappa1
%           Ten run with 'BL' and 'Info' and 'ksi']
VARModelStr = 'BL';  %'CHOL';  %    Baseline (BL) or standard high-frequency identification (CHOL)
VARModelShock = 'Info';  %'MonPol'; Monetary policy shock 'MonPol' or information shock 'Info' 
VARModelInfoShock = 'ksi'; %or 'ksi' or 'g'; 'betta'; %Which model shock is used for info shock

switch_OptType = 'SA';  %Simulated annealing

%To test how sensitive is the monetary policy shock to financial frictions
flagKappa0  = 0;            %Sets kappa to 0; 1 if kappa=0, 0 if kappa kappa>0
if flagKappa0==1 
    switch_OptType = '';    %overwrite: no optimization, uses the starting values (starting values need to set to optimum) 
end;

flagRun = 0;        %1 to run impulse response matching
flagPlot = 1;       %1 to plot impulse response functions

%Structure of variable names collected from the model (dynare output)
varsCell       =   {'y1_ann';  'N'; 'Y'; 'P'; 'ebp_ann'}; %
varsCellNames  =   {'1-year rate'; 'S&P500'; 'GDP'; 'GDP deflator'; 'EBP'}; %

%Structure of variable names collected from the VAR (irf output)
varsCell(:,2)         =   {'gs1';  'logsp500'; 'loggdp'; 'loggdpdef'; 'ebp'}; %
varsCellNames(:,2)  =   {'i'; 'Q'; 'Y';'P';  'ebp'}; %
rangePlot = [-0.2 -3  -0.3  -0.2   -0.1      %maximal range use on plots
              0.4  2   0.2   0.2   0.1];
%The variable names in the model and in the VAR should match
          
switch VARModelShock
    case 'Info'
        switch VARModelInfoShock
            case 'ksi'
                varsExoCell      =   {'e_ksi'};  %
                varsExoCellNames   =   {'e_{\xi}'}; %
            case 'g'
                varsExoCell      =   {'e_g'};  %
                varsExoCellNames   =   {'e_{g}'}; %
            case 'betta'
                varsExoCell      =   {'e_betta'};  %
                varsExoCellNames   =   {'e_{\beta}'}; %
            case 'a'
                varsExoCell      =   {'e_a'};  %
                varsExoCellNames   =   {'e_{a}'}; %
        end;                
    otherwise
        varsExoCell      =   {'e_ir'};  %
        varsExoCellNames   =   {'e_{i}'}; %
end;

switch VARModelShock
    case 'Info'
        varsExoCell(:,2)      =   {'e_info'};  %
        varsExoCellNames(:,2)   =   {'e_{info}'}; %
    otherwise
        varsExoCell(:,2)      =   {'e_ir'};  %
        varsExoCellNames(:,2)   =   {'e_{i}'}; %
end;

NVars  =   size(varsCell,1);
NVarsExo = size(varsExoCell,1);

%Loading impulse responses
load Irfs/us1_ff4_hf_sp500_hf_1984m02-2016m12_chol_mom12.mat
CHOLMean = mom.mean./1e2;
CHOLVar  = mom.var./1e4;
load Irfs/us1_ff4_hf_sp500_hf_1984m02-2016m12_sgnm2_mom12.mat
BLMean = mom(1,1).mean./1e2;
BLVar = mom(1,1).var./1e4;
BLInfoMean = mom(1,2).mean./1e2;
BLInfoVar = mom(1,2).var./1e4;

IrfLen=length(BLMean)/5;

switch VARModelStr
    case 'CHOL'
        irf2Mean = CHOLMean;
        irf2Var  = CHOLVar;
    case 'BL'
        switch varsExoCell{:,2}
            case 'e_ir'
                irf2Mean = BLMean;
                irf2Var  = BLVar;                
            case 'e_info'
                irf2Mean = BLInfoMean;
                irf2Var  = BLInfoVar;                
        end;
end;
irf2MeanMatr = reshape(irf2Mean(:),[IrfLen 5]);
irf2 = irf2MeanMatr;
irf2Low2Sd = irf2 - 2*reshape(sqrt(diag(irf2Var)),[IrfLen 5]);
irf2High2Sd = irf2 + 2*reshape(sqrt(diag(irf2Var)),[IrfLen 5]);
irf2Low1Sd = irf2 - reshape(sqrt(diag(irf2Var)),[IrfLen 5]);
irf2High1Sd = irf2 + reshape(sqrt(diag(irf2Var)),[IrfLen 5]);

if flagRun==1
%Create VAR irf matrix
%kk=2;
%run CreateIrfMatr.m;

%Run the full dynare code first
cd dynare;
options_.order = 1;
options_.irf = IrfLen;
options_.periods = 0;        
dynare 'GKHH_dynare' noclearall;

%Then run a code to optimize the parameters
%choose variables, set minimum and maximum values and transform variables
%for optimization
switch VARModelShock
    case 'Info'
        switch VARModelInfoShock
            case 'ksi'
                ParamsOptCell = {'sigma_ksi','rhoKsi'};
                XXMin         = [1e-4       1e-5     ];
                XXMax            = [0.01       1-1e-5];
                TransType = {'exp','tan'};
            case 'g'
                ParamsOptCell = {'sigma_g','rhoG'};
                XXMin         = [1e-4       1e-5     ];
                XXMax            = [0.01       1-1e-5];
                TransType = {'exp','tan'};
            case 'betta'
                ParamsOptCell = {'sigma_betta','rhoShockBetta'};
                XXMin         = [1e-4       1e-5     ];
                XXMax            = [0.01       1-1e-5];
                TransType = {'exp','tan'};
            case 'a'
                ParamsOptCell = {'sigma_a','rhoA'};
                XXMin         = [1e-4       1e-5     ];
                XXMax            = [0.01       1-1e-5];
                TransType = {'exp','tan'};
        end;
    otherwise
        ParamsOptCell = {'sigma_ir','rhoShockIr','kappa','gam_P','gam'};
        XXMin         = [5e-4          0.5      1e-3       1e-5     1e-5];
        XXMax         = [0.001      1-1e-5     1        1-1e-5    1-1e-5];
        TransType = {'exp','tan','exp','tan','tan'};
end;

%Set initial values for optimization
switch VARModelStr
    case 'CHOL'
        XX0           = [0.00085     0.65      0.005     0.95     0.9];   %CHOL
    case 'BL'
        switch VARModelShock
            case 'Info'
                switch VARModelInfoShock
                    case 'ksi'
                        XX0           = [0.00055     0.85];  %BLInfo
                    case 'g'
                        XX0           = [0.0045     0.7582];  %BLInfo
                    case 'betta'
                        XX0           = [0.0008     0.6240];  %BLInfo
                    case 'a'
                        XX0           = [0.00055     0.85];  %BLInfo
                end;
            otherwise
                if flagKappa0==0
                    XX0           = [0.0008    0.6371    0.0452    0.0000    0.8649];   
                else  %setting kappa to a low number (use the optimized values) 
                    XX0           = [0.0008    0.6371    1e-3    0.0000    0.8649];   
                end;
        end;
end;
NParams = length(ParamsOptCell);

%Create structure for parameters
ObjFuncParams = v2struct(varsCell,varsCellNames,varsExoCell,varsExoCellNames,VARModelStr,...
            NVars,NVarsExo,ParamsOptCell,NParams,oo_,switch_OptType,TransType,... 
            VARModelShock,VARModelInfoShock,rangePlot,flagKappa0,...
            irf2Mean,irf2Var,irf2MeanMatr,IrfLen,irf2,irf2Low2Sd,irf2High2Sd,irf2Low1Sd,irf2High1Sd);

switch switch_OptType
    case 'SA'  %simulated annealing
        optionsSA = saoptimset('Display','iter','InitialTemperature',[1 1],'TemperatureFcn',@temperatureexp,'StallIterLimit',250,'HybridFcn',@fmincon,'TolFun',1e-8);
        YY0 = var_transform(XX0,TransType,'input'); 
        YYMin = var_transform(XXMin,TransType,'input'); 
        YYMax = var_transform(XXMax,TransType,'input'); 
        [YY,objfunc,exitf] = simulannealbnd(@(ZZ) ObjFunc(ZZ,ObjFuncParams),YY0,YYMin,YYMax,optionsSA);
        XX = var_transform(YY,TransType,'output');
    otherwise  %no optimization
        XX = XX0;
end;

%Display the results
disp(ParamsOptCell);
disp(XX);

% Creates model impulse responses
run ../CreateModelIrf.m
ObjFuncParams.irf1 = irf1;
cd ..;

%saving results
eval(['save data/ObjFuncParams_' VARModelStr '_' VARModelShock '_' VARModelInfoShock '_flagKappa0' num2str(flagKappa0) '.mat ObjFuncParams;']);
end; %end of optimization

%plotting (assumes optimization was run at least once)
if flagPlot==1
%plotting results
    VARModelStrCell = {     'BL',       'BL',   'CHOL'};
    VARModelShockCell = {   'MonPol',   'Info', 'MonPol'};
    VARModelInfoShockCell={ 'ksi',      'ksi',  'ksi'};
    NNFigures = length(VARModelStrCell);

    pathStr    =   'figures/';  
    extStr    = '';  %any extra string at the end of the figure file name

    for kk=1:NNFigures
    %loading results
        eval(['load data/ObjFuncParams_' VARModelStrCell{kk} '_' VARModelShockCell{kk} '_' VARModelInfoShockCell{kk} '_flagKappa00.mat ObjFuncParams;']);
%plotting results
        plot_figures(ObjFuncParams,kk,'VAR','b--','ff');
        plot_figures(ObjFuncParams,kk,'Model','k','gg');
        if strcmp(VARModelStrCell{kk},'BL') &&  strcmp(VARModelShockCell{kk},'MonPol')
            eval(['load data/ObjFuncParams_' VARModelStrCell{kk} '_' VARModelShockCell{kk} '_' VARModelInfoShockCell{kk} '_flagKappa01.mat ObjFuncParams;']);
            plot_figures(ObjFuncParams,kk,'Model','r.','hh');
            ff=get(gca,'Children');
            ff=flipud(ff);
            legend([ff(3) ff(4) ff(5)],'VAR','Model','Model,\kappa=0');
        else
            ff=get(gca,'Children');
            ff=flipud(ff);
            legend([ff(3) ff(4)],'VAR','Model');
        end;        
        plot_figures(ObjFuncParams,kk,'ZeroLine','k','qq');
    end;

    for kk=1:NNFigures
         eval(['figures_save(pathStr,extStr,''GKHH_VAR_' VARModelStrCell{kk} '_' VARModelShockCell{kk} '_' VARModelInfoShockCell{kk} ''',kk);']);
    end;
end;