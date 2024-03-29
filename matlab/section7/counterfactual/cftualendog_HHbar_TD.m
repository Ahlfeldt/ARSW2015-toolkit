%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program is new, but it builds on the replication directory     %%%
%%% It perfoms illustrative straightforward counterfactuals to          %%%
%%% highlight the potential of the model                                %%%
%%% It generates various maps and graphs for didactic purposes          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We use the smodendog.m programme to find the values of the target
%%% variables wages, floor space prices, and commercial floor space shares,
%%% based on which other endogenous outcomes are created
%%% These counterfactual correspond to the closed-city case

%%% The illustrative counterfactuals include a complete car ban, and a
%%% productivity boost to East Berlin. Other counterfactuals can be
%%% perfomed in the same way, by changing other primitives

%%% We compare these counterfactuals under endogenous agglomeration forces 
%%% to the counterfactuals with exogenous fundamentals

%%%%%%%%%%%%%%%
%%% Car ban %%%
%%%%%%%%%%%%%%%

% Clear all data 
clear all;
clc;
clf;
colormap default;
format bank;
close all;

% Load data set generated by the quantification in cftualprep_end_TD.m
load('data/output/endogcftual_prep_big_TD');

% Load counterfactual transport cost matrix
load('data/input/ttpublic_2006_ren');    

%%% Solve for conterfactual equilibrium without car travel %%%%%%%%%%%%%%%%
%%% Takes about 20 minutes to run on 2024 LSE-RLABTS15 Server %%%%%%%%%%%%%
tic;
[Cendog_cft,Cucprob_cft,HH_cft,Cconverge_cft,Ubar06_cft]=smodendog(param,fund,ttpub06,nobs06);
Cwage_cft=Cendog_cft(:,1); Cvv_cft=Cendog_cft(:,2); Ctheta_cft=Cendog_cft(:,3);
CYY_cft=Cendog_cft(:,4); CQ_cft=Cendog_cft(:,5); Cq_cft=Cendog_cft(:,6);
CHM_cft=Cendog_cft(:,7); CHR_cft=Cendog_cft(:,8);
wtime=toc;
display('>>>> Time to solve for equilibrium (mins) <<<<');
wtime./60
% Save data for illustration
save('data/output/endogcftual_TT_HHbar','Cendog_cft','Cucprob_cft','HH_cft','Cconverge_cft','Ubar06_cft','');

%%% Store results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute percentage changes 
UpccTT_end_HHbar = (Ubar06_cft./Ubar06-1).*100                              % Expected utility falls by 2.75%, vs. 5.7% with exogenous fundamentals!
% Compute percentage changes and map
pccTT_end_HHbar = (Cendog_cft ./ Cendog -1).*100;                           % Indeed, total employment virtually constant
pccTT_end_HHbar(pccTT_end_HHbar>100)=100;   % Truncate changes at -50% and +100% for a clarity of the map
pccTT_end_HHbar(pccTT_end_HHbar<-50)=-50;   % Truncate changes at -50% and +100% for a clarity of the map



%%% Inspect results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Residence employment
display('>>>> Mapping change in residence employment <<<<');
pccTT_end_HHbar_CHR=pccTT_end_HHbar(:,8);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_HHbar_CHR,'Car ban, closed city: % change in residence employment','figs/endcft','MAP_pcc_Transport_HHbar_HR_HHbar') 
% Remaining residents leave suburbs and move towards the center 

% Workplace employment
display('>>>> Mapping change in workplace employment <<<<');
pccTT_end_HHbar_CHM=pccTT_end_HHbar(:,7);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_HHbar_CHM,'Car ban, closed city: % change in workplace employment','figs/endcft','MAP_pcc_Transport_HM_HHbar') 
% Firms are croweded out of the center and become relatively more present
% in suburbs - blocks become less specialized

% floorspace prices
display('>>>> Mapping change in floor space prices <<<<');
pccTT_end_HHbar_CQ=pccTT_end_HHbar(:,5);
pccTT_end_HHbar_Cq=pccTT_end_HHbar(:,6);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_HHbar_CQ,'Car ban, closed city: % change in floor space prices','figs/endcft','MAP_pcc_Transport_Q_HHbar') 
% Center becomes more expensive

%%% Comparison to counterfactual with exogenous fundamentals %%%%%%%%%%%%%%

% Residence employment
load('data/output/TTcft_exog','UpccTT_exog','pccTT_exog')  
pccTT_exog_CHR=pccTT_exog(:,8);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(emprsd06>0),pccTT_end_HHbar_CHR(emprsd06>0));
hold on;
scatter(distCBDr(emprsd06>0),pccTT_exog_CHR(emprsd06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Change in residential employment');
title('Car ban, closed city: % change in residential employment');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
%saveas(gcf,'figs/endcft/FIG_TT_OpenClosed_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_TT_EndExog_HR.png','-dpng','-r300');

% Compare change in workplace employment gradients %%%%%%%%%%%%%%%%%%%%%%%%
pccTT_exog_CHM=pccTT_exog(:,7);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(empwpl06>0),pccTT_end_HHbar_CHM(empwpl06>0));
hold on;
scatter(distCBDr(empwpl06>0),pccTT_exog_CHM(empwpl06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Change in workplace employment');
title('Car ban, closed city: % change in workplace employment');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
% saveas(gcf,'figs/endvsexog/FIG_TT_EndExog_HM.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_TT_EndExog_HM.png','-dpng','-r300');

% Compare change floor space gradients %%%%%%%%%%%%%%%%%%%%%%%%
pccTT_exog_CQ=pccTT_exog(:,5);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(emprsd06>0),pccTT_end_HHbar_CQ(emprsd06>0));
hold on;
scatter(distCBDr(emprsd06>0),pccTT_exog_CQ(emprsd06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Change in residential floor space price');
title('Car ban, closed city: % change in residential floor space price');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
% saveas(gcf,'figs/endvsexog/FIG_TT_EndExog_HM.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_TT_EndExog_Q.png','-dpng','-r300');

% Compare change floor space gradients %%%%%%%%%%%%%%%%%%%%%%%%
pccTT_exog_Cq=pccTT_exog(:,5);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(empwpl06>0),pccTT_end_HHbar_Cq(empwpl06>0));
hold on;
scatter(distCBDr(empwpl06>0),pccTT_exog_Cq(empwpl06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Change in commercial floor space price');
title('Car ban, closed city: % change in commercial floor space price');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
% saveas(gcf,'figs/endvsexog/FIG_TT_EndExog_HM.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_TT_EndExog_q.png','-dpng','-r300');

% Compare change in wage gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pccTT_exog_Cw=pccTT_exog(:,1);
pccTT_end_HHbar_Cw=pccTT_end_HHbar(:,1);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(empwpl06>0),pccTT_end_HHbar_Cw(empwpl06>0));
hold on;
scatter(distCBDr(empwpl06>0),pccTT_exog_Cw(empwpl06>0));
hold off;
xlabel('Distance from the CBD');
ylabel('Change in wage');
title('Car ban, closed city: % change in wage');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
% saveas(gcf,'figs/endvsexog/FIG_TT_EndExog_HM.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_TT_EndExog_w.png','-dpng','-r300');

% Save counterfactual results for later comparison
save('data/output/TTcft_end_HHbar','UpccTT_end_HHbar','pccTT_end_HHbar'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Productivity increase in East Berlin %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear all data 
clear all;
clc;
clf;
colormap default;
format bank;
close all;

% Load data set generated by the quantification in cftualprep_end_TD.m
load('data/output/endogcftual_prep_big_TD');

% Adjust productivity in the East %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dummyeast = dummywestr.*-1 +1;
aCtf = a06.*(1+0.1.*dummyeast);
fund(:,1)=aCtf;

%%% Solve for counterfactual equilibrium and save endogenous outcomes %%%%%
tic;
[Cendog_cft,Cucprob_cft,HH_cft,Cconverge_cft,Ubar06_cft]=smodendog(param,fund,tt06,nobs06);
wtime=toc;
display('>>>> Time to solve for equilibrium (mins) <<<<');
wtime./60

% Save data for illustration
save('data/output/endogcftual_AE_HHbar','Cendog_cft','Cucprob_cft','HH_cft','Cconverge_cft','Ubar06_cft');

%%% Store results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute percentage changes and map
pccAE_end_HHbar = (Cendog_cft ./ Cendog -1).*100;                           % Compute percentage changes
UpccAE_end_HHbar = (Ubar06_cft./Ubar06-1).*100                              % Expected utility increases 3.1% (compared to 3.8% under exogenous fundamentals)
pccAE_end_HHbar(pccAE_end_HHbar>100)=100;   % Truncate changes at -50% and +100% for a clarity of the map
pccAE_end_HHbar(pccAE_end_HHbar<-50)=-50;   % Truncate changes at -50% and +100% for a clarity of the map

% Save counterfactual results for later comparison
save('data/output/AEcft_end_HHbar','UpccAE_end_HHbar','pccAE_end_HHbar'); 

%%% Inpect results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wages
display('>>>> Mapping change in wage <<<<');
pccAE_end_HHbar_Cw=pccAE_end_HHbar(:,1);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_HHbar_Cw,'East Berlin productivity gain: % change in wage','figs/endcft','MAP_pcc_EastProd_w_HHbar') 
% 

% Residence employment
display('>>>> Mapping change in residence employment <<<<');
pccAE_end_HHbar_CHR=pccAE_end_HHbar(:,8);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_HHbar_CHR,'East Berlin productivity, closed city: % change in residence employment','figs/endcft','MAP_pcc_EastProd_HHbar_HR_HHbar') 
% Remaining residents leave suburbs and move towards the center 

% Workplace employment
display('>>>> Mapping change in workplace employment <<<<');
pccAE_end_HHbar_CHM=pccAE_end_HHbar(:,7);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_HHbar_CHM,'East Berlin productivity: % change in workplace employment','figs/endcft','MAP_pcc_EastProd_HM_HHbar') 
% Firms are croweded out of the center and become relatively more present
% in suburbs - blocks become less specialized

% floorspace prices
display('>>>> Mapping change in floor space prices <<<<');
pccAE_end_HHbar_CQ=pccAE_end_HHbar(:,5);
pccAE_end_HHbar_Cq=pccAE_end_HHbar(:,6);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_HHbar_CQ,'East Berlin productivity: % change in floor space prices','figs/endcft','MAP_pcc_EastProd_Q_HHbar') 
% Center becomes more expensive

%%% Comparison to counterfactual with exogenous fundamentals %%%%%%%%%%%%%%

% Load counterfactual effects under exogenous fundamentals
load('data/output/AEcft_exog','pccAE_exog');   

% Compare change in wage gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pccAE_exog_Cw=pccAE_exog(:,1);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(empwpl06>0),pccAE_end_HHbar_Cw(empwpl06>0));
hold on;
scatter(distWALLrunning(empwpl06>0),pccAE_exog_Cw(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in wage');
title('East Berlin productivitiy upgrade: % change in wage');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_EastProd_EndExog_w.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_EastProd_EndExog_w.png','-dpng','-r300');
% Wages respond more to change in fundamental productivity due to
% agglomeration effects. More positively in east, more negatively in the
% west. But notice the increase in wages with endogenous agglomeration
% forces close to the border due to agglomeration spillovers from the east.

% Compare change in workplace employment gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pccAE_exog_CHM=pccAE_exog(:,7);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(empwpl06>0),pccAE_end_HHbar_CHM(empwpl06>0));
hold on;
scatter(distWALLrunning(empwpl06>0),pccAE_exog_CHM(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in wage');
title('East Berlin productivitiy upgrade: % change in workplace employment');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_EastProd_EndExog_HM.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_EastProd_EndExog_HM.png','-dpng','-r300');
% Workplace employment responds more to change in fundamental productivity due to
% agglomeration effects

% Compare change in residence employment gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pccAE_exog_CHR=pccAE_exog(:,8);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(emprsd06>0),pccAE_end_HHbar_CHR(emprsd06>0));
hold on;
scatter(distWALLrunning(emprsd06>0),pccAE_exog_CHR(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in residence employment');
title('East Berlin productivitiy upgrade: % change in residence employment employment');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_EastProd_EndExog_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_EastProd_EndExog_HR.png','-dpng','-r300');
% Residence employment responds more to change in fundamental productivity due to
% agglomeration effects

% Compare change in residential floor space prices %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pccAE_exog_CQ=pccAE_exog(:,5);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(emprsd06>0),pccAE_end_HHbar_CQ(emprsd06>0));
hold on;
scatter(distWALLrunning(emprsd06>0),pccAE_exog_CQ(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in floor space prices');
title('East Berlin productivitiy upgrade: % change in floor space prices');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_EastProd_EndExog_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_EastProd_EndExog_HQ.png','-dpng','-r300');
% Price changes follow a smooth trend across the border due to spillovers.
% With endogenous agglomeration effects the price responses are generally
% larger due to the aplication effect on productivity coming from
% employment relocation.

% Compare change in commercial floor space prices %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pccAE_exog_Cq=pccAE_exog(:,6);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(empwpl06>0),pccAE_end_HHbar_Cq(empwpl06>0));
hold on;
scatter(distWALLrunning(empwpl06>0),pccAE_exog_Cq(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in commercial floor space prices');
title('East Berlin productivitiy upgrade: % change in commercial floor space prices');
legend('Endogenous agglomeration forces','Exogenous fundamentals','Location','best');
%saveas(gcf,'figs/endvsexog/FIG_EastProd_EndExog_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/endvsexog/FIG_EastProd_EndExog_Hq.png','-dpng','-r300');
% Price changes follow a smooth trend across the border due to spillovers.
% With endogenous agglomeration effects the price responses are generally
% larger due to the aplication effect on productivity coming from
% employment relocation.

display('<<< Counterfactuals with endogenous agglomeration force and exogenous total employment successfully completed >>>')


