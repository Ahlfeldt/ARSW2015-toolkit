%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB programme file for the toolkit for           %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This programme is based on a programme in the replication directory %%%
%%% This version has been commented                                     %%%
%%% This version focuses on 1986 data to reduce memory use              %%%
%%% This version refers to Epsi (instead of theta) as the Frechet       %%%
%%% shape parameter (consistent with notations in the paper             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration of the parameter epsilon for Section 6 of the paper;
clear all;
clc;
clf;
colormap default;
format shortG;
close all;

% Declaring (creating) objects as globals so they can be accessed and modified by any function
global alpha beta epsilon kappaeps kappa lambda delta eta rho eps;
global lambda0 delta0 eta0 rho0 ThetaD ThetaP ThetaR ThetaA Theta0 UBD UBP UBR LBD LBP LBR;
global ftD ftP ftR ftA ND LD WD LP WP LR WR fwestd fwestr bzk86dw bzk86rw
global obsvar36 tt36 nobs36 omega_36 omega36 wage36 a36 b36;
global obsvar06 tt06 nobs06 omega_06 omega06 wage06 a06 b06;
global obsvar86dw tt86dw nobs86dw omega_86dw omega86dw wage86dw a86dw b86dw;
global obsvar86rw tt86rw nobs86rw omega_86rw omega86rw wage86rw a86rw b86rw;
global IICBDdw IICBDrw distCBDdw distCBDrw IIBZKdw IIBZKrw;
global Xdwmat Ydwmat Xrwmat Yrwmat Xdwcut Ydwcut Xrwcut Yrwcut;
global varlwdata varlwmod lbwdata lbwmod;

% *******************;
% **** Load Data ****;
% *******************;

load('data/input/prepdata_big_TD86');

% ********************************;
% **** Read Bezirke wage data ****;
% ********************************;

[bzkwge]=csvread('data/input/wageworker1986.csv');                          % Reads data from csv containin wages by Bezirke
bzkwge=bzkwge(:,2);                                                         % Keep column 2 containing numerical wages
lbwdata=log(bzkwge);                                                        % Generate log wages
lbwdata=lbwdata-mean(lbwdata);                                              % Demean log wages
varlwdata=var(lbwdata);                                                     % Compute variance of log wages, our empirical moment

% ***********************;
% **** Random Number ****;
% ***********************;

% Set default random number stream;
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% Tolerance;
eps = 1e-6;                                                                 % Define tolerance level for the numerical procedure

% *************************;
% **** Initializations ****;
% *************************;

% Alpha and beta;
alpha=0.80;                                                                 % Set parameter values to values from the literature
beta=0.75;                                                                  % Set parameter values to values from the literature
kappaeps=0.07;                                                              % Set commuting decay to reduced-form estimate

obsvar86rw=zeros(nobs86rw,4);                                                % Generate a vector of key variables used by transformed wages solver
obsvar86rw(:,1)=floor86rw; obsvar86rw(:,2)=empwpl86rw; obsvar86rw(:,3)=emprsd86rw; 

omega_86rw=zeros(nobs86rw,1);                                               % Generate starting values for transformed wages to speed up convergence
omega_86rw(empwpl86rw>0)=(((1-alpha)./floor86rw(empwpl86rw>0)).^((1-alpha)./alpha)).*alpha;

% *************************;
% **** CALIBRATE OMEGA ****;
% *************************;

% We use our comegaopt function to calibrate for 1986 reunification;
[omega86rw,cprob86rw,wconverge86rw,HMC86rw,wgap86rw] = comegaoptO(obsvar86rw,tt86rw,nobs86rw,omega_86rw);
omega_86rw=omega86rw;
% We take a moment to enjoy what we have achieved by inspecting transformed wages on a map
lomega86rw = log(omega86rw+1);
RESULT = MAPIT('../shapefile/WestBerlin4matlab',lomega86rw,'Transformed wages','figs','calibration/MAP_omega_86') 

display('>>>> Convergence wage 86rw Details <<<<');
display([wconverge86rw wgap86rw]);
display('>>>> Time to solve for wages (mins) <<<<');

display('***********************************************');
display('**** REUNIFICATION ONE-STEP GMM ESTIMATION ****');
display('***********************************************');

% COMMUTING ESTIMATION;
LBD = 2;                        
UBD = 24;
Epsi0 = 4;
% FMINCON;
% options = optimset('display','iter','TolFun',eps,'TolX',eps,'Algorithm','interior-point');
% [ThetaD,ll,EXITFLAG,OUTPUT,LAMBDA] = ...
%   fmincon('cdensityoptren',ThetaD0,[],[],[],[],LBD,UBD,[],options);
% PATTERN SEARCH;
options = psoptimset('Display','iter','TolFun',eps,'TolX',eps);
[Epsi,ll,exitflag,output] = patternsearch(@cdensityoptren,Epsi0,[],[],[],[],LBD,UBD,[],options); 
    % Patternsearch is part of the Global Optimization Toolbox. It is a
    % flexible search algorithm that minimizes a value of an objective
    % function
        % We define this function in cdensityoptren.m, it will compute the
        % residual sum of squares between the model-implied log wage
        % variance across Bezirke and the log variance observed in data
    % Patternsearches for the parameter value for Epsi, taking Epsi0 as
        % a starting value

display('>>>> Estimated Parameter <<<<');
display('>>>> [epsilon] <<<<');
epsilon

Ewage86rw=wage86rw(wage86rw>0);                                             % Model-consistent 1986 wages for blocks with positive wages
varlwage86rw=var(log(Ewage86rw));                                           % Compute block-level variance on log wages

display('>>>> Variance Log Bezirke Wages 1986 <<<<');
varlwdata

display('>>>> Variance Log Block Wages 1986 <<<<');
varlwage86rw

% We take a moment to enjoy what we have achieved by inspecting adjusted wages on a map
RESULT = MAPIT('../shapefile/WestBerlin4matlab',wage86rw,'Adjusted wages','figs','calibration/MAP_adjWage_86') 

% **********************************;
% **** Assumed Parameter Values ****;
% **********************************;

epsilon=round(epsilon.*100)./100;
kappa=kappaeps./epsilon;

% *******************;
% **** SAVE DATA ****;
% *******************;

keep varlwage86rw;

save('data/output/blockvarlwage86rw_big');

display('>>>> One-step estimation of epsilon completed successfully <<<<');
