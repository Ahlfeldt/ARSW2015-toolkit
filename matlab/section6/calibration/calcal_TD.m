%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This version of this program is part of the replication directory   %%%
%%% This version has been commented                                     %%%
%%% This version has marginally modified for greater clarity            %%%
%%% This version generates various maps for didactic purposes           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This script executes the seuqential procedure for the quantification 
%%% with exogenous fundamentals

clear all;
clc;
clf;
colormap default;
format shortG;
close all;

% Declaring (creating) objects as globals so they can be accessed and modified by any function
global alpha beta kappa epsilon eps;
global fwestd fwestr;
global  A06  ;
global  B06  ;
global  tt06  ;
global  empwpl06  ;
global  emprsd06  ;
global  nobs06  ;

% *******************;
% **** Load Data ****;
% *******************;

load('data/input/prepdata_big_TD');

% ***********************;
% **** Random Number ****;
% ***********************;

% Set default random number stream;
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% Tolerance;
eps = 1e-6;

% ********************;
% **** PARAMETERS ****;
% ********************;

% Alpha and beta;
alpha=0.80;                                                                 % Set parameter values to values from the literature 
beta=0.75;                                                                  % Set parameter values to values from the literature
% Kappa and epsilon;
kappaeps=0.07;                                                              % Reduced-form estimate of commuting decay
epsilon=6.83;                                                               % Preference heterogeneity estimated in optimepsilon_TD86.m
kappa=kappaeps./epsilon;                                                    % Implied iceberg commuting cost parameter
% Gamma function;
epsfrac=(epsilon-1)./epsilon;                                               % gamma from Eq. (9)
gamfun=gamma(epsfrac);                                                      % gamma from Eq. (9)

% *************************;
% **** Initializations ****;
% *************************;

obsvar06=zeros(nobs06,4);                                                   % Generate a vector of key variables, floor space prices, workplace employment, residence employment, geographic area
obsvar06(:,1)=floor06; obsvar06(:,2)=empwpl06; obsvar06(:,3)=emprsd06; obsvar06(:,4)=area06;

wage_06=zeros(nobs06,1);                                                    % Bring wage guesses into plausible range to speed up convergence 
wage_06(empwpl06>0)=(((1-alpha)./floor06(empwpl06>0)).^((1-alpha)./alpha)).*alpha;
    % This is equation (12) sovlved for wages after setting A = 1
    % This brings guessed values into a plauisble range and should help
    % with the speed of the solver
    
A_06=zeros(nobs06,1); A_06(empwpl06>0)=1;                                   % Starting values for production fundamentals set to 1
B_06=zeros(nobs06,1); B_06(emprsd06>0)=1;                                   % Starting values for residential fundamentals set to 1

% Modern Bezirke;
[modbzk06] = modbezirk(bzk06,nobs06);                                       % Simple function that applies a cross-walk to generate identifier of 12 modern Bezirke from 23 historic Bezirke

% *************************;
% **** CALIBRATE WAGES ****;
% *************************;

tic;                                                                        % Start timer
% Calibrate wages and productivities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We use a solver similar to the one in comegaoptO.m to recover adjusted
% wages using workplace employment, residence employment, and bilateral travel times in  (Eq. S.44) in the supplement
% Given transformed wages, this version of the solver also recovers
% adjusted poductivity using Eq. (27)
[wage06,A06,cprob06,wconverge06,HMC06,wgap06] = comegaoptC(obsvar06,tt06,nobs06,wage_06);
display('>>>> Convergence wage 06 Details <<<<');
display([wconverge06 wgap06]);
wtime=toc;
display('>>>> Time to solve for wages (mins) <<<<');
wtime./60

% Check if we corectly recover workplace employment
scatter(empwpl06,HMC06);                                                    % Yes, slope = 1
xlabel('Workplace employment in data');
ylabel('Workplace employment in model');
[b, stats] = regress(HMC06, empwpl06) % 

% We take a moment to enjoy what we have achieved by inspecting adjusted
display('>>>> Map production fundamentals <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',A06,'Adjusted productivity','figs','calibration/MAP_A_06') 
    % Observe large production fundamentals in central area near the Mitte
    % CBD and the City West. These are the locations with highest floor
    % space prices
display('>>>> Map floor space prices <<<<');
lfloor06 = log(floor06); 
RESULT = MAPIT('../shapefile/Berlin4matlab',lfloor06,'Log floor space prices','figs','calibration/MAP_lfloorspaceprices_06') 
    % These are also locations with high wages
display('>>>> Map adjusted wages <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',wage06,'Adjusted wages','figs','calibration/MAP_wage_06') 

display('>>>> Map workplace employment densities <<<<');
lEmpWplDens06 = log( empwpl06./area06 + 1);
RESULT = MAPIT('../shapefile/Berlin4matlab',lEmpWplDens06,'Log workplace employment density in data','figs','calibration/MAP_lEmpWplDens_06') 
display('>>>> Map residential employment densities <<<<');
lEmpRsdDens06 = log( emprsd06./area06 + 1);
RESULT = MAPIT('../shapefile/Berlin4matlab',lEmpRsdDens06,'Log residence employment density in data','figs','calibration/MAP_lEmpRsdDens06_06') 

% *****************************;
% **** CALIBRATE AMENITIES ****;
% *****************************;

tic;
% Calibrate amenities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We use a simple solver that compute commuting market access and then
% solves for adjusted amenities using Eq. S.47.
[B06,CMA06,HRS06] = camen(obsvar06,tt06,nobs06,wage06,cprob06);
btime=toc;
display('>>>> Time to solve for amenities (mins) <<<<');
btime./60

% Check if we corectly recover residence employment
H_obs = sum(emprsd06)                                                      % Total number of workers, needed for comparion since HRS06 are block shares while emprsd06 are block counts
scatter(emprsd06,HRS06.*H_obs)                                             % Yes, slope = 1
xlabel('Residence employment in data');
ylabel('Residence employment in model');
[b, stats] = regress(HRS06.*H_obs, emprsd06)                               % Yes, slope = 1
% Notice that we have residence employment SHARE, not levels.
% This is because we are not adjusting the level to ensure that we have the
% right number of workers in our open city. In other words, our amenities
% are identified up to a multiplicative scale factor

% We take a moment to enjoy what we have achieved by inspecting adjusted
% amenities on a map
display('>>>> Map adjusted amenities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',B06,'Adjusted amenity','figs','calibration/MAP_B_06') 
display('>>>> Map commuter market access <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',CMA06,'Commuter market access','figs','calibration/MAP_CMA06_06') 


% *******************************************************;
% **** ADJUST LEVELS OF PRODUCTIVITIES AND AMENITIES ****;
% *******************************************************;

%%% If we want to use the recovered primitives in the counterfactual
%%% analysis, we need to make sure that A and B are in the right levels. To
%%% this end, we need to rescale wages to ensure that A have a geometric
%%% mean of one and rescale B so that the population in model (Phi) matches
%%% the population in the data. Then, A and B will correspond to the values
%%% recovered by cmdexog.m. 
%%% The function below makes this adjustment. It will change the outputs
%%% compared to the original replacation directory. You do not need to use
%%% the function if you are not interested in the levels of A and B and
%%% your do not intend to use A and B from the sequential procedure as
%%% input into the smodexpg.m solver

% Update productivities and amenities
[A06,B06,wage06] = calcal_adj_TD(obsvar06,tt06,nobs06,A06,B06);

% *****************************;
% **** TOTAL WORKER INCOME ****;
% *****************************;

% Total worker income for 2006;
% We use a solver that computes expected worker income at residence as wage 
% at workplace weighted by conditional commuting probabilities 
[vv06] = expincome(obsvar06,tt06,nobs06,wage06,B06);

% We take a moment to enjoy what we have achieved by inspecting adjusted
display('>>>> Map expected total income <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',vv06,'Expected income','figs','calibration/MAP_vv06_06') 


% ******************************************;
% **** CALIBRATE DENSITY OF DEVELOPMENT ****;
% ******************************************;

% Calibrate density of development for 2006;
[V06,L06,theta06] = cdensity(obsvar06,A06,wage06,vv06,nobs06);              % This is a simple function that uses Eqs. S29-S31 to solve for density of development, total floor space, and share of residential floor space

% We take a moment to enjoy what we have achieved by inspecting adjusted
display('>>>> Map density of development <<<<');
lV06 = log(V06+1);
RESULT = MAPIT('../shapefile/Berlin4matlab',lV06,'Density of development','figs','calibration/MAP_lV06_06') 
display('>>>> Map total floor space <<<<');
lL06 = log(L06+1);
RESULT = MAPIT('../shapefile/Berlin4matlab',lL06,'Total floor space','figs','calibration/MAP_lL06_06') 
display('>>>> Map commercial floor space share <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',theta06,'Share of commercial floor space','figs','calibration/MAP_theta06_06') 

V06rw=V06(fwestr);
L06rw=L06(fwestr);


% ******************;
% **** Clean Up ****;
% ******************;

clear Xd Xdmat Xdw Xdwmat Xr Xrmat Xrw Xrwmat Yd Ydmat Ydw Ydwmat Yr Yrmat Yrw Yrwmat;
clear distvec06 distvec36;
clear dist06 dist36 dist86dw dist86rw;
clear distvec06 distvec36;
clear Xd Xdmat Xdw Xdwmat Xr Xrmat;
clear Xrw Xrwmat Yd Ydmat Ydw Ydwmat Yr Yrmat Yrw Yrwmat;
clear data06 data36 data86d data86dw data86r data86rw;

% *******************;
% **** Save Data ****;
% *******************;
clear cprob06
LD06 = L06;                                                                 % Replicating so that the objet (floor space stock) will be found by the equilibrium solver
save('data/output/calcal_big');
% end;

display('>>>> Sequential quantification succesfully completed <<<<');
