%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A variant of this program is part of the replication directory      %%%
%%% This version has been commented                                     %%%
%%% This version has modified for greater clarity and reduced memory    %%%
%%% use; it focuses on objects that are relevant for the 2006           %%%
%%% quantification.                                                     %%%
%%% This version has been extended to illustrate the breakdown into     %%%    
%%% fundamental and endogenous productivities and amenities             %%%
%%% productivities and amenities recovered in Section 6                 %%% 
%%% This version has been extended to provide a comparison to the       %%%
%%% productivities and amenities recovered in Section 6                 %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script prepares the data for counterfactuals;
% This entails quantifying the model

% Inputs 
    % cmodexog.m 
    % cdensityE.m
    % cprod.m
    % cres.m
    % ubar.m
    
    
% Clear inputs   
clear all;
clc;
clf;
colormap default;
format shortG;
close all;

% Declaring (creating) objects as globals so they can be accessed and modified by any function
global alpha beta kappa epsilon lambda delta rho eta eps;
global fwestd fwestr;
global  A06 ;
global  B06 ;
global  tt06 ;
global  empwpl06 ;
global  emprsd06 ;
global  nobs06 ;
global IICBDdw IICBDrw;
global WP NP LP cutoff;

% *******************;
% **** Load Data ****;
% *******************;

%if bigdata==1;
load('data/input/prepdata_big_TD');                                               % This is the basic data prior to quantification containing objects relevant to this teaching directory
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
alpha=0.80;                                                                 % Set input share of labour in production 
beta=0.75;                                                                  % Set expenditure share on floor space    

% ****************************************;
% *** LOAD POOLED PARAMETER ESTIMATES ****;
% ****************************************;

%if bigdata==1;
load('data/input/roptimis_all_big');                             % Load estimated structural parameters
%end;

% Extract parameter values from Theta vector from the GMM estimation of the
% full model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappaeps=EThetaA(1);                                                        % Commuting decay (nu)
epsilon=EThetaA(2);                                                         % Preference heterogeneity
kappa=kappaeps./epsilon;                                                    % Commuting cost     
lambda=EThetaA(3);                                                          % Density elasticity of productivity
delta=EThetaA(4);                                                           % Density elasticity of residential amenity
eta=EThetaA(5);                                                             % Productivity decay    
rho=EThetaA(6);                                                             % Residenital decay        

% Use parameter values from Section 6 instead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If desired, uncomment the following lines to use the same parameter values
% as in Section 6 with exogenous fundamentals. The code will then recover
% the same productivities and amenities A and B as in Section 6, with the
% only difference being that they are broken down into an endogenous and an
% exogenous component.
%{
% Kappa and epsilon;
kappaeps=0.07;                                                              % Reduced-form estimate of commuting decay
epsilon=6.83;                                                               % Preference heterogeneity estimated in optimepsilon_TD86.m
kappa=kappaeps./epsilon;                                                    % Implied iceberg commuting cost parameter
%}

display('>>>> Parameters <<<<');
display('[kappaeps epsilon lambda delta eta rho]');
display([kappaeps epsilon lambda delta eta rho]);

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

% Create placeholders for fundamental productivities and amenities, will be
% updated by cmodexog.m programme
A_06=zeros(nobs06,1); A_06(empwpl06>0)=1;
B_06=zeros(nobs06,1); B_06(emprsd06>0)=1;
    % This is equation (12) sovlved for wages after setting A = 1
    % This brings guessed values into a plauisble range and should help
    % with the speed of the solver
    
% ************************************************;
% **** CALIBRATE PRODUCTIVITIES AND AMENITIES ****;
% ************************************************;

tic;
%%% Calibrate productivities and amenities for 2006; %%%%%%%%%%%%%%%%%%%%%%
% Also save a range of endogenous outcomes
[A06,B06,wage06,ucprob06,vv06,HMC06,HRC06,CMA06,Ephi06,HH06,ABconverge06,mAgap06,mBgap06] = cmodexog(obsvar06,tt06,nobs06,A_06,B_06);
display('[Converge Agap Bgap]');
[ABconverge06 mAgap06 mBgap06]

HRS06=emprsd06./sum(emprsd06);                                              % Compute residence employment shares
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

% *********************************;
% **** Production Fundamentals ****;
% *********************************;

%%% Break down productivities in endogenous and fundamental component %%%%%
[a06,Ups06] = cprod(obsvar06,tt06,nobs06,A06);

% **********************************;
% **** Residential Fundamentals ****;
% **********************************;

%%% Break down amenities in endogenous and fundamental component %%%%%
[b06,Ome06] = cres(obsvar06,tt06,nobs06,B06);

% ******************;
% **** Clean Up ****;
% ******************;

clear wage_36 wage_86dw wage_06 wage_86rw;
clear A_36 A_86dw A_06 A_86rw;
clear B_36 B_86dw B_06 B_06rw;
clear Xd Xdmat Xdw Xdwmat Xr Xrmat Xrw Xrwmat Yd Ydmat Ydw Ydwmat Yr Yrmat Yrw Yrwmat;
clear distvec06 distvec36;
clear dist06 dist36 dist86dw dist86rw;
clear distvec06 distvec36;
clear Xrw Xrwmat Yd Ydmat Ydw Ydwmat Yr Yrmat Yrw Yrwmat;
clear data06 data36 data86d data86dw data86r data86rw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialisations for counterfactual %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Read parameter values into vector as input into ussmodendog.m %%%%%%%%%
param=[kappa epsilon lambda delta eta rho];

%%% Read fundamentals and endogenous variables into vector as solver input 
fund=zeros(nobs06,10);
fund(:,1)=a06; fund(:,2)=b06; fund(:,3)=V06; fund(:,4)=area06;
fund(:,5)=floor06; fund(:,6)=empwpl06; fund(:,7)=emprsd06; 
fund(:,8)=wage06; fund(:,9)=vv06; fund(:,10)=theta06;

%%% Compute reservation utility level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In this open-city application, we keep utility constant at the
%%% reservation utility level and let total employment adjust 
[Ubar06] = ubar(B06,param,fund,tt06,nobs06)                                 % New function based on code in replication directory

%%% Solve for initial equilibrium and save endogenous outcomes %%%%%%%%%%%%
[Cendog,Cucprob,HH,CUbar,Cconverge]=ussmodendog(param,fund,tt06,nobs06,Ubar06);
Cwage=Cendog(:,1); Cvv=Cendog(:,2); Ctheta=Cendog(:,3);
CYY=Cendog(:,4); CQ=Cendog(:,5); Cq=Cendog(:,6);
CHM=Cendog(:,7); CHR=Cendog(:,8);


% *******************;
% **** Save Data ****;
% *******************;

%if bigdata==1;
save('data/output/endogcftual_prep_big_TD');
%end;

 display('>>>> Data preparation for counterfactuals completed <<<<');

% ********************;
% **** Inspection ****;
% ********************;

%%% Inspect productivities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('>>>> Map endogenous adjusted productivities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',Ups06,'Endogenous productivity','figs/endcft','MAP_Ups06') 
% There is a obvious core-periphery pattern with endogenous productivities
% being highest in the centers with the greatest employment densities
% (Mitte, Potsdamer Platz, Kurfuerstendamm)

display('>>>> Map fundamental adjusted productivities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',a06,'Exogenous productivity','figs/endcft','MAP_a06')
% They appear spatially random, model has explained the core-periphery
% pattern via endogenous density effect

% Scatter plot of gradients
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(A06>0),A06(A06>0));
hold on;
scatter(distCBDr(a06>0),a06(a06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Productivity');
title('Productivity gradient');
legend('Total','Exogenous','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_TotalExog_Productivities.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endcft/FIG_TotalExog_Productivities.png','-dpng','-r300');
% Scatter plot confirms that the peak productivities are explained by
% endogenous productivity
mdl = fitlm(distCBDr(A06>0),log(A06(A06>0)))                                % Total productivity declines by 1.7% per km from CBD
mdl = fitlm(distCBDr(a06>0),log(a06(a06>0)))                                % Exogenous productivity increases by 0.2% per km from CBD

%%% Inspect amenities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('>>>> Map endogenous adjusted amenities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',Ome06,'Endogenous amenity','figs/endcft','MAP_Ome06') 
% We see a doughnut of endogenous amenities in the areas just outside the
% commercial areas and within (and just outside) the S-Bahn ring where
% residential densities are high

display('>>>> Map exogenous adjusted amenities <<<<');
RESULT = MAPIT('../shapefile/Berlin4matlab',b06,'Exogenous amenity','figs/endcft','MAP_b06') 
% Fundamental amenities are mostly falt in the central parts of teh city,
% but there are often high values in suburban locations close to green
% spaces and water bodies


% Scatter plot of gradients
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(B06>0),B06(B06>0));
hold on;
scatter(distCBDr(b06>0),b06(b06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Amenity');
title('Amenity gradient');
legend('Total','Exogenous','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_TotalExog_Amenities.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endcft/FIG_TotalExog_Amenities.png','-dpng','-r300');
% Scatter plot confirms that the peak amenities clsoe to the CBD are 
% explained by endogenous amenity. The exogenous is upward sloping,
% consistent with natural amenities in the suburbs
mdl = fitlm(distCBDr(B06>0),log(B06(B06>0)))                                % Total amenity declines by 0.4% per km from CBD
mdl = fitlm(distCBDr(b06>0),log(b06(b06>0)))                                % Exogenous amenity increases by 1.8% per km from CBD

% ************************************************************;
% **** Comparson to exogenous fundamentals from Section 6 ****;
% ************************************************************;

% Notice that we use exactly the same procedure as in Section 6. The only
% difference is that we now use different parameter values from the GMM
% esimtation of the model with endogenous productivities and amenities. The
% differences between the recovered primitives from Section 6 and those
% recovered here for Section 7 are, thus, directly informative of
% sensitivity to the parameter values. The most notable difference in
% parameter values is that nu=epsilon*kappa, at 0.0987 (vs. 0.07) is
% larger. Epsilon, at 6.6941 (vs. 6.83), is marginally smaller.

Acft = A06; Bcft = B06; LDcft  = LD06 ;                                     % Rename fundamentals
load('data/output/exogcftual_prep_big_TD','A06','B06', 'LD06');                    % Load fundamentals from sequential procedure
% Compare productivities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
scatter(A06,Acft); 
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Section 6 estimates');
ylabel('Section 7 estimates');
title('Productivities')
mdl = fitlm(Acft,A06)
saveas(gcf,'figs/endvsexog/FIG_Sec6vsSec7_Productivities.png');

% Compare amenities
clf;
scatter(B06,Bcft);
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Section 6 estimates');
ylabel('Section 7 estimates');
title('Amenities')
mdl = fitlm(Bcft,B06)
saveas(gcf,'figs/endvsexog/FIG_Sec6vsSec7_Amenities.png');
% Under the larger kappa*epsilon estimate (which implies greater commuting
% costs), we rationalize the data with greater amenities

% Compare floor space stock
clf;
scatter(LD06,LDcft);
hold on;
xlims = xlim;
plot(xlims, xlims);
hold off;
xlabel('Sequential procedure');
ylabel('Simultaneous procedure');
title('Floor space stock')
mdl = fitlm(LDcft,LD06)
saveas(gcf,'figs/endvsexog/FIG_Sec6vsSec7_FloorSpace.png');

display('>>>> Comparison of cases with and without endogenous agglomeration forces completed <<<<');
