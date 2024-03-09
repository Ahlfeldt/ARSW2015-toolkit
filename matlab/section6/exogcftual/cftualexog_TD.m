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
%%% It has been modified illustrate the convergence path and test for   %%% 
%%% the correct recovery of endogenous outcomes                         %%%
%%% It perfoms some illustrative straightforward counterfactuals to     %%%
%%% highlight the potential of the model                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We use the smodexog.m programme to find the values of the targe
%%% variables wages, floor space prices, and commercial floor space shares,
%%% based on which other endogenous outcomes are created
%%% We show that we correctly recover the initial values of the endogenous
%%% variables
%%% The illustrative counterfactuals include a complete car ban, and a
%%% productivity boost to East Berlin. Other counterfactuals can be
%%% perfomed in the way, by changing other primitives

% Clear all data 
clear all;
clc;
clf;
colormap default;
format bank;
close all;

% Load data set generated by the quantification in calcal_TD.m 
% This is the quantification from the sequential procedure. It is not the
% quantification used in the paper. It will only work if you exectute the
% added calcal_adj_TD.m program from calcal_TD.m as this will ensure that
% the productivities and amenities are in the levels required to match the
% intial population. 
% load('data/output/calcal_big');

% Load data set generated by the quantification in cftualprep_TD.m
% Notice that this is quantification used in the paper
load('data/output/exogcftual_prep_big_TD')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What happens if we ban cars in the entire city?                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data/input/ttpublic_2006_ren');                                             % Load public transport travel time matrix

% Now run the solver with the new transport matrix %%%%%%%%%%%%%%%%%%%%%%%%
[CendogCft,CucprobCft,HHCft,CconvergeCft,cpathCft,UbarCft]=smodexog(param,fund,ttpub06,nobs06);
% Save data for quick loading 
save('data/output/exogcftual_prep_big_TD_TT','CendogCft','CucprobCft','HHCft','CconvergeCft','cpathCft','UbarCft'); 

% Let's inspect convergence
close all
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
plot(cpath(:,7),(cpath(:,1)),'b-','LineWidth',2) 
hold on;
plot(cpath(:,7),(cpath(:,2)),'r--','LineWidth',2) 
hold on;
plot(cpath(:,7),(cpath(:,4)),'color',[0,0.5,0],'LineStyle',':','LineWidth',2) 
legend({'Wages','Floor space prices','Land use'},'Location','best')
xlabel('Iterations') 
ylabel('Max. log gap')
set(fig,'Papersize',[6 3]);
print(fig,'figs/exogcft/FIG_cftualexog_maxLD_TT.png','-dpng','-r300');

% Change in expected utility
UpccTT = (UbarCft ./ Ubar - 1)   .* 100                                    % Utility falls by 5.7%
% Compute percentage changes 
pccTT = (CendogCft ./ Cendog -1).*100; 
pccTT(pccTT>100)=100;   % Truncate changes at -50% and +100% for a clarity of the map
pccTT(pccTT<-50)=-50;   % Truncate changes at -50% and +100% for a clarity of the map

% Residence employment
display('>>>> Mapping change in residence employment <<<<');
pccTT_CHR=pccTT(:,8);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_CHR,'Car ban: % change in residence employment','figs/exogcft','MAP_pcc_Transport_HR') 
% Remaining residents leave suburbs and move towards the center 

% Workplace employment
display('>>>> Mapping change in workplace employment <<<<');
pccTT_CHM=pccTT(:,7);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_CHM,'Car ban: % change in workplace employment','figs/exogcft','MAP_pcc_Transport_HM') 
% Firms are croweded out of the center and become relatively more present
% in suburbs - blocks become less specialized

% #floorspace prices
display('>>>> Mapping change in floor space prices <<<<');
pccTT_CQ=pccTT(:,5);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_CQ,'Car ban: % change in floor space prices','figs/exogcft','MAP_pcc_Transport_Q') 
% Center becomes more expensive

% Save counterfactual results for later comparison
pccTT_exog = pccTT;
UpccTT_exog = UpccTT;
save('data/output/TTcft_exog','UpccTT_exog','pccTT_exog');   

display('<<< Transport conterfactual completed >>>')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What happens if we make East Berlin more attractive?               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear all data 
clear all;
clc;
clf;
colormap default;
format bank;
close all;

% Load data set generated by the quantification in calcal.m 
load('data/output/exogcftual_prep_big_TD')

% Adjust productivity in the East %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dummyeast = dummywestr.*-1 +1;
ACtf = A06.*(1+0.1.*dummyeast);
fund(:,1)=ACtf;
% Use solver to solve for the counterfactual equilibrium %%%%%%%%%%%%%%%%%%
[CendogCft,CucprobCft,HHCft,CconvergeCft,cpathCft,UbarCft]=smodexog(param,fund,tt06,nobs06);
% Save data for quick loading 
save('data/output/exogcftual_prep_big_TD_AE','CendogCft','CucprobCft','HHCft','CconvergeCft','cpathCft','UbarCft') 

% Let's inspect convergence
close all
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
plot(cpath(:,7),(cpath(:,1)),'b-','LineWidth',2) 
hold on;
plot(cpath(:,7),(cpath(:,2)),'r--','LineWidth',2) 
hold on;
plot(cpath(:,7),(cpath(:,4)),'color',[0,0.5,0],'LineStyle',':','LineWidth',2) 
legend({'Wages','Floor space prices','Land use'},'Location','best')
xlabel('Iterations') 
ylabel('Max. log gap')
set(fig,'Papersize',[6 3]);
print(fig,'figs/exogcft/FIG_cftualexog_maxLD_TT.png','-dpng','-r300');

% Change in expected utility
UpccAE_exog = (UbarCft ./ Ubar - 1)   .* 100                                % Utility increaes by 3.8%                                                   
% Compute percentage changes and map
pccAE = (CendogCft ./ Cendog -1).*100; 
pccAE(pccAE>100)=100;   % Truncate changes at -50% and +100% for a clarity of the map
pccAE(pccAE<-50)=-50;   % Truncate changes at -50% and +100% for a clarity of the map

% Workplace employment
display('>>>> Mapping change in workplace employment <<<<');
pccAE_CHM=pccAE(:,7);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_CHM,'Eastern renewal: % change in workplace employment','figs/exogcft','MAP_pcc_EastRenewal_HM') 

% Residence employment
display('>>>> Mapping change in residence employment <<<<');
pccAE_CHR=pccAE(:,8);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_CHR,'Eastern renewal: % change in residence employment','figs/exogcft','MAP_pcc_EastRenewal_HR') 

% floorspace prices
display('>>>> Mapping change in floor space prices <<<<');
pccAE_CQ=pccAE(:,5);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_CQ,'Eastern renewal: % change in floor space prices','figs/exogcft','MAP_pcc_EastRenewal_Q') 

% wages
display('>>>> Mapping change in wages <<<<');
pccAE_Cwage=pccAE(:,1);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_Cwage,'Eastern renewal: % change in wages','figs/exogcft','MAP_pcc_EastRenewal_w') 

% Save counterfactual results for later comparison
pccAE_exog = pccAE;
save('data/output/AEcft_exog','pccAE_exog','UpccAE_exog');   

display('<<< Productivity conterfactual completed >>>')

display('>>>> Counterfactuals with exogenous fundamentals completed <<<<');