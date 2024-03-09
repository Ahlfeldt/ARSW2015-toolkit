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
%%% Several descriptive excercises have been added for didactic         %%%
%%% purposes                                                            %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This script prepares the data for counterfactuals;
%%% This entails quantifying the model and solving for the initial
%%% equilibrium
%%% This script also provides a comparison to the fundamentals inverted
%%% using the sequential procedure. It confirms that, upon rescaling of the
%%% fundamentals from the sequential procedure, the fundamentals are
%%% identical

% This script uses following functions 
    % cmodexog.m 
    % cdensityE.m
    % smodexog.m
  
    
% Clear inputs    
clear all;
clc;
clf;
colormap default;
format shortG;
close all;

% Declaring (creating) objects as globals so they can be accessed and modified by any function
global alpha beta kappa epsilon eps;
global fwestd fwestr;
global A06 
global B06 
global tt06 
global empwpl06 
global emprsd06 
global nobs06 

% *******************;
% **** Load Data ****;
% *******************;

load('data/input/prepdata_big_TD');                                               % This is the basic data prior to quantification
%end;

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

wage_06=zeros(nobs06,1);                                                     % Bring wage guesses into plausible range to speed up convergence 
wage_06(empwpl06>0)=(((1-alpha)./floor06(empwpl06>0)).^((1-alpha)./alpha)).*alpha;
    % This is equation (12) sovlved for wages after setting A = 1
    % This brings guessed values into a plauisble range and should help
    % with the speed of the solver
    
% Create placeholders for fundamental productivities and amenities, will be
% updated by cmodexog.m programme
A_06=zeros(nobs06,1); A_06(empwpl06>0)=1;                                   % Productivities
B_06=zeros(nobs06,1); B_06(emprsd06>0)=1;                                   % Amenities

% ************************************************;
% **** CALIBRATE PRODUCTIVITIES AND AMENITIES ****;
% ************************************************;

tic;
%%% Calibrate productivities and amenities for 2006; %%%%%%%%%%%%%%%%%%%%%%
% Also save a range of endogenous outcomes
[A06,B06,wage06,ucprob06,vv06,HMC06,HRC06,CMA06,Ephi06,HH06,ABconverge06,mAgap06,mBgap06] = cmodexog(obsvar06,tt06,nobs06,A_06,B_06);
display('[Converge Agap Bgap]');
[ABconverge06 mAgap06 mBgap06]

HRS06=emprsd06./sum(emprsd06);
%HRS86rw=emprsd86rw./sum(emprsd86rw);
btime=toc;
display('>>>> Time to solve for amenities (mins) <<<<');
btime./60

% *******************************;
% **** CALIBRATE LAND MARKET ****;
% *******************************;

%%%% Calibrate density of development for 2006; %%%%%%%%%%%%%%%%%%%%%%%%%%%
[V06,LD06,LM06,LR06,theta06] = cdensityE(obsvar06,A06,wage06,vv06,nobs06);

% *****************;
% **** OUTPUT *****;
% *****************;

%%% Output for 2006; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y06=A06.*(empwpl06.^alpha).*(LM06.^(1-alpha));                              % We use the production function in Eq. (10)

% ********************;
% **** Inspection ****;
% ********************;

% We take a moment to enjoy what we have achieved by inspecting adjusted
display('>>>> Map adjusted productivities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',A06,'Adjusted productivity','figs/exogcft','MAP_A06') 
display('>>>> Map adjusted amenities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',B06,'Adjusted amenity','figs/exogcft','MAP_B06') 
display('>>>> Map commuter market access <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',CMA06,'Commuter market access','figs/exogcft','MAP_CMA06') 
display('>>>> Map adjusted wage <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',wage06,'Adjusted wage','figs/exogcft','MAP_wage06') 


% ******************;
% **** Clean Up ****;
% ******************;

clear wage_36 wage_86dw wage_06 wage_86rw;
clear A_36 A_86dw A_06 A_86rw;
clear B_36 B_86dw B_06 B_06rw;
clear Xd Xdmat Xdw Xdwmat Xr Xrmat Xrw Xrwmat Yd Ydmat Ydw Ydwmat Yr Yrmat Yrw Yrwmat;
clear distvec06 distvec36;
clear distvec06 distvec36;
clear dist36 dist86dw dist06 dist86rw;
clear Xrw Xrwmat Yd Ydmat Ydw Ydwmat Yr Yrmat Yrw Yrwmat;
clear data06 data36 data86d data86dw data86r data86rw;

% ****************;
% **** Inputs ****;
% ****************;

% Prepare inputs for equilibrium solvers
    % Notice htat only A05 (productivities) and B06 (amenities) are
    % fundamentals. The other variables are endogenous outcomes that are solved
    % within the model. We do not need those values to solve for the
    % equilibrium. However, using them speeds up convergence. 
    % We can replace some (or all) of the endogenous variables with onevec
    % (a dummy of ones) and observe how it takes longer for convergence to
    % be reached. Still, we recover values that are close to those observed
    % in the data, confirming that inversion and solver work hand-in-glove.
onevec = ones(nobs06,1);
param=[epsilon kappaeps];
fund=zeros(nobs06,12);
fund(:,1)=A06; fund(:,2)=B06; 
fund(:,3)= V06; % onevec; % V06; % Here, we replace the solved equilibrium values with vectors of ones to make it more demanding for the solver to find the equilibrium
fund(:,4)=area06;
fund(:,5)=  floor06; % onevec; % floor06; 
fund(:,6)= empwpl06; % onevec; % empwpl06; 
fund(:,7)=  emprsd06; % onevec; % emprsd06;
% Commercial, residential and total floor space are endogenous objects solved within the model. Commercial and residential floor space are being generated by applying the commercial floor space share to total floor space
fund(:,8)=  LD06.*theta06;  % onevec.*0.5; % L06.*theta06; 
fund(:,9)= LD06.*(1-theta06); % onevec.*0.5; % L06.*(1-theta06); 
fund(:,10)=LD06; %  % Total floor space is inverted and treated as exogenous   
fund(:,11)= wage06;  % onevec; % wage06; 
fund(:,12)=  vv06; % onevec; % vv06;

% ***********************************************;
% **** Run solver and save baseline outcomes ****;
% ***********************************************;

% Use solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Cendog,Cucprob,HH,Cconverge,cpath,Ubar]=smodexog(param,fund,tt06,nobs06);
wtime=toc;
display('>>>> Time to solve for equilibrium (mins) <<<<');
wtime./60

% Save the endogenous variables generated by smodexog.m
Cwage=Cendog(:,1); Cvv=Cendog(:,2); Ctheta=Cendog(:,3);
CYY=Cendog(:,4); CQ=Cendog(:,5); Cq=Cendog(:,6);
CHM=Cendog(:,7); CHR=Cendog(:,8); Crent=Cendog(:,9);

% Update wage input for solver for faster convergence
fund(:,5)= Crent;  % onevec; % wage06; 
fund(:,11)= Cwage;  % onevec; % wage06; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inspect outputs: Do we recover workplace employment correctly?        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do we recover workplace employment corretly?
clf;
scatter(CHM ,empwpl06) 
mdl = fitlm(CHM,empwpl06) %  'CHM ~ empwpl06-1'
% Yes!

% Do we recover residence amenities correctly?
scatter(CHR ,emprsd06) 
mdl = fitlm(CHR,emprsd06) %  'CHM ~ empwpl06-1'
% Yes!

diary on;
display('>>>> Check 2006 Converged <<<<');
Cconverge
display('>>>> Check 2006 Reproduce Initial Equilibrium <<<<');
display('Total employment');
display('[sum(empwpl) sum(CHM) sum(emprsd) sum(CHR) sum(Y) sum(CY)]');
[sum(empwpl06) sum(CHM) sum(emprsd06) sum(CHR)  sum(CYY)] % sum(Y06)
display('[empwpl; emprsd; floorQ; floorq; Y; wage]');
[max(abs(empwpl06-CHM)); max(abs(emprsd06-CHR)); max(abs(floor06-CQ)); max(abs(floor06-Cq)); max(abs(wage06-Cwage))] % ; max(abs(Y06-CYY))
diary off;

% *******************;
% **** Save Data ****;
% *******************;

save('data/output/exogcftual_prep_big_TD');
%end;

display('>>>> Data for counterfactuals saved <<<<');

% *******************************************;
% **** Comparson to sequential procedure ****;
% *******************************************;

Acft = A06; Bcft = B06; LDcft  = LD06 ;                                     % Rename fundamentals
clear A06 B06 LD06;
load('data/output/calcal_big','A06','B06', 'L06');                        % Load fundamentals from sequential procedure
% Compare productivities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
scatter(A06,Acft); 
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Sequential procedure');
ylabel('Simultaneous procedure');
title('Productivities')
mdl = fitlm(Acft,A06)
% We obtain exactly the same productivities as with the sequential
% procedure. Notice that this will not be the case if you do not call the 
% calcal_adj_TD.m from calcal_TD.m  

% Compare amenities
clf;
scatter(B06,Bcft);
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Sequential procedure');
ylabel('Simultaneous procedure');
title('Amenities')
mdl = fitlm(Bcft,B06)
% We obtain exactly the same amenities as with the sequential
% procedure. Notice that this will not be the case if you do not call the 
% calcal_adj_TD.m from calcal_TD.m  

% Compare floor space stock
clf;
scatter(L06,LDcft);
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Sequential procedure');
ylabel('Simultaneous procedure');
title('Floor space stock')
mdl = fitlm(LDcft,L06)
% We obtain exactly the same floor space stock as with the sequential
% procedure. Notice that this will not be the case if you do not call the 
% calcal_adj_TD.m from calcal_TD.m  

%%% Fundamentals strongly correlated, but they scale differently. Need to
%%% use the fundamentals inverted by cmodexog.m in counterfactuals

lB06 = log(B06+1); lBcft = log(Bcft+1);
scatter(lB06,lBcft);
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Sequential procedure');
ylabel('Simultaneous procedure');
title('Ln amenities')
mdl = fitlm(lBcft,lB06)

%%% Let's inspect the geometric means of productivities
EA06 = A06(A06>0);
geomean(EA06)           % Sequential procedure
EAcft = Acft(Acft>0);
geomean(EAcft)          % Simultaneous procedure
%%% Sequential procedure generates higher productivity levels since we
%%% normalize wages to their geometric means and not productivities

%%% Let's inspect the geometric means of amenities
EB06 = B06(B06>0);
geomean(EB06)           % Sequential procedure
EBcft = Bcft(Bcft>0);
geomean(EBcft)          % Simultaneous procedure
%%% Sequential proceure generates amenities normalized to a geometric mean
%%% of one. We require a higher level of amenities to generate a Phi (total
%%% employment in model) that matches total population.

display('>>>> Data preparation for counterfactuals with exogenous fundamentals completed <<<<');

