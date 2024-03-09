%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The orginal version of this file was part of the replication        %%%
%%% directory to Ahlfeldt, Redding, Sturm, Wolf (2015)                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This version has been commented                                     %%%
%%% This version has marginally modified for greater clarity            %%%
%%% by Gabriel M. Ahlfeldt in 01/2024                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration of the parameter epsilon for Section 6 of the paper;
clear all;
clc;
clf;
colormap default;
format shortG;
close all;

% User==1 is Steve Windows desktop;
% User==2 is Steve Mac Pro desktop;
% User==5 is Berlin Humboldt server;

% bigdata==1 is the full dataset;

% Select user and data;
user=1;
bigdata=1;
    
%if user==1;
 %   cd 'E:\Data\Gabriel\Ahfleldt\_QSE\ARSW\FinalBerlin\matlab\section6\optimepsilon';
%elseif user==2;
%    cd '/Users/reddings/Dropbox/SRDS/FinalBerlin/matlab/section6/optimepsilon';
%elseif user==5;
%    cd '/work/eodwiwi/FinalBerlin/matlab/section6/optimepsilon';
%end;

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

if bigdata==1;
load('data/prepdata_big_TD');
end;

% ********************************;
% **** Read Bezirke wage data ****;
% ********************************;

[bzkwge]=csvread('data/wageworker1986.csv');                                % Reads data from csv containing wages by Bezirke
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

obsvar06=zeros(nobs06,4);                                                   % Generate a vector of key variables used by transformed wages solver
obsvar06(:,1)=floor06; obsvar06(:,2)=empwpl06; obsvar06(:,3)=emprsd06; obsvar06(:,4)=area06;
    % now we introduce a vector of variables that are central for solving
    % for epsilon
        % normalized floor space prices
        % employment at workplace
        % employment at residence
        % block area

%obsvar86rw=zeros(nobs86rw,4);                                                % Generate a vector of key variables used by transformed wages solver
%obsvar86rw(:,1)=floor86rw; obsvar86rw(:,2)=empwpl86rw; obsvar86rw(:,3)=emprsd86rw; obsvar86rw(:,4)=area86rw;

omega_06=zeros(nobs06,1);                                                    % Generate starting values for transformed wages to speed up convergence
omega_06(empwpl06>0)=(((1-alpha)./floor06(empwpl06>0)).^((1-alpha)./alpha)).*alpha; 
    % This is equation (12) sovlved for wages after setting A = 1
    % This brings guessed values into a plauisble range and should help
    % with the speed of the solver

%omega_86rw=zeros(nobs86rw,1);                                               % Generate starting values for transformed wages to speed up convergence
%omega_86rw(empwpl86rw>0)=(((1-alpha)./floor86rw(empwpl86rw>0)).^((1-alpha)./alpha)).*alpha;

%wage_06=zeros(nobs06,1);
%wage_06(empwpl06>0)=(((1-alpha)./floor06(empwpl06>0)).^((1-alpha)./alpha)).*alpha;

%wage_86rw=zeros(nobs86rw,1);
%wage_86rw(empwpl86rw>0)=(((1-alpha)./floor86rw(empwpl86rw>0)).^((1-alpha)./alpha)).*alpha;

% *************************;
% **** CALIBRATE OMEGA ****;
% *************************;

% We use our comegaopt function to calibrate trasnformed wages for 2006;
[omega06,cprob06,wconverge06,HMC06,wgap06] = comegaoptO(obsvar06,tt06,nobs06,omega_06);
omega_06=omega06;
display('>>>> Convergence wage 06 Details <<<<');
display([wconverge06 wgap06]);

% We take a moment to enjoy what we have achived by inspecting transformed wages on a map
lomega_06 = log(omega_06+1);
RESULT = MAPIT('..\shapefile\Berlin4matlab',lomega_06,'Transformed wages','figs','MAP_omega') 

% We use our comegaopt function to calibrate for 1986 reunification;
%[omega86rw,cprob86rw,wconverge86rw,HMC86rw,wgap86rw] = comegaopt(obsvar86rw,tt86rw,nobs86rw,omega_86rw);
%omega_86rw=omega86rw;
% We take a moment to enjoy what we have achieved by inspecting transformed wages on a map
%lomega86rw = log(omega86rw+1);
%RESULT = MAPIT('..\shapefile\WestBerlin4matlab',lomega86rw,'Transformed wages','figs','MAP_omega_86') 

%display('>>>> Convergence wage 86rw Details <<<<');
%display([wconverge86rw wgap86rw]);
%display('>>>> Time to solve for wages (mins) <<<<');
%{
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
RESULT = MAPIT('..\shapefile\WestBerlin4matlab',wage86rw,'Adjusted wages','figs','MAP_adjWage_86') 

% **********************************;
% **** Assumed Parameter Values ****;
% **********************************;

epsilon=round(epsilon.*100)./100;
kappa=kappaeps./epsilon;

% *******************;
% **** SAVE DATA ****;
% *******************;

keep varlwage86rw;

if bigdata==1;
save('../../data/blockvarlwage86rw_big');
end;
%}
display('>>>> File Completed Successfully <<<<');
