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
%%% It has been modified illustrate the convergence path and test for   %%% 
%%% the correct recovery of endogenous outcomes                         %%%
%%% It perfoms illustrative straightforward counterfactuals to          %%%
%%% highlight the potential of the model                                %%%
%%% It generates various maps and graphs for didactic purposes          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We use the ussmodendog.m programme to find the values of the target
%%% variables wages, floor space prices, and commercial floor space shares,
%%% based on which other endogenous outcomes are created
%%% These counterfactual correspond to the open-city case

%%% The illustrative counterfactuals include a complete car ban, and a
%%% productivity boost to East Berlin. Other counterfactuals can be
%%% perfomed in the same way, by changing other primitives

%%% We compare counterfactuals to the closed-city case

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
[Cendog_cft,Cucprob_cft,HH_cft,CUbar_cft,Cconverge,cpath]=ussmodendog(param,fund,ttpub06,nobs06,Ubar06);
wtime=toc;
display('>>>> Time to solve for equilibrium (mins) <<<<');
wtime./60
% Save data for illustration
save('data/output/endogcftual_TT','Cendog_cft','Cucprob_cft','HH_cft','CUbar_cft','Cconverge','cpath');

%%% Store results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change in total employment
delta_HH_TT_end = HH_cft - HH                                               % city loses 240k workers (280k workers with exogenous fundamentals)
% Compute percentage changes and map
pccTT_end = (Cendog_cft ./ Cendog -1).*100; 
pccTT_end(pccTT_end>100)=100;   % Truncate changes at -50% and +100% for a clarity of the map
pccTT_end(pccTT_end<-50)=-50;   % Truncate changes at -50% and +100% for a clarity of the map


%%% Inspect convergence path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
plot(cpath(:,7),(cpath(:,1)),'b-','LineWidth',2) 
hold on;
plot(cpath(:,7),(cpath(:,2)),'r--','LineWidth',2) 
hold on;
plot(cpath(:,7),(cpath(:,4)),'color',[0,0.5,0],'LineStyle',':','LineWidth',2) 
legend({'Wages','Floor space prices','Land use'},'Location','best')
xlabel('Iterations') 
ylabel('Max. log gap')
saveas(gcf,'figs/endcft/FIG_cftualTT_maxLD.png');

close all
plot(cpath(:,7),(cpath(:,6)),'color','blue','LineStyle','-','LineWidth',2) 
legend({'MSLE'},'Location','best')
xlabel('Iterations') 
ylabel('Log MSLE')
saveas(gcf,'figs/endcft/FIG_cftualTT_MSLE.png');

%%% Inspect results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Residence employment
display('>>>> Mapping change in residence employment <<<<');
pccTT_end_CHR=pccTT_end(:,8);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_CHR,'Car ban, open city: % change in residence employment','figs/endcft','MAP_pcc_Transport_HR') 
% Remaining residents leave suburbs and move towards the center 

% Workplace employment
display('>>>> Mapping change in workplace employment <<<<');
pccTT_end_CHM=pccTT_end(:,7);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_CHM,'Car ban, open city: % change in workplace employment','figs/endcft','MAP_pcc_Transport_HM') 
% Firms are croweded out of the center and become relatively more present
% in suburbs - blocks become less specialized

% floorspace prices
display('>>>> Mapping change in floor space prices <<<<');
pccTT_end_CQ=pccTT_end(:,5);
pccTT_end_Cq=pccTT_end(:,6);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_CQ,'Car ban, open city: % change in floor space prices','figs/endcft','MAP_pcc_Transport_Q') 
% Center becomes more expensive

% Wage
display('>>>> Mapping change in wage <<<<');
pccTT_end_Cw=pccTT_end(:,1);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccTT_end_Cw,'Car ban, open city: % change in wage','figs/endcft','MAP_pcc_Transport_w') 
% Most wages fall, except for some areas in teh suburbs where employment
% denisty increases as the remaining residents move to more central areas

% Save counterfactual results for later comparison
pccTT_end = pccTT_end;
save('data/output/TTcft_end','delta_HH_TT_end','pccTT_end');   


%%% Comparison to closed-city case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Residence employment
load('data/output/TTcft_end_HHbar','UpccTT_end_HHbar','pccTT_end_HHbar'); 
pccTT_end_HHbar_CHR=pccTT_end_HHbar(:,8);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(emprsd06>0),pccTT_end_CHR(emprsd06>0));
hold on;
scatter(distCBDr(emprsd06>0),pccTT_end_HHbar_CHR(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the CBD');
ylabel('Change in residential employment');
title('Car ban: % change in residential employment');
legend('Open city','Closed city','Location','best');
%saveas(gcf,'figs/endcft/FIG_TT_OpenClosed_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_TT_OpenClosed_HR.png','-dpng','-r300');
% In the open-city case the suburbs loose more residence employment, and
% there are fewer people left in the city, so the center gains less

% Workplace employment
pccTT_end_HHbar_CHM=pccTT_end_HHbar(:,7);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(empwpl06>0),pccTT_end_CHM(empwpl06>0));
hold on;
scatter(distCBDr(empwpl06>0),pccTT_end_HHbar_CHM(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the CBD');
ylabel('Change in workplace employment');
title('Car ban: % change in workplace employment');
legend('Open city','Closed city','Location','best');
%saveas(gcf,'figs/endcft/FIG_TT_OpenClosed_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_TT_OpenClosed_HM.png','-dpng','-r300');
% In the open-city case there are fewer people left in the city, so the 
% center gains less

% Wage
pccTT_end_HHbar_Cw=pccTT_end_HHbar(:,1);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(empwpl06>0),pccTT_end_Cw(empwpl06>0));
hold on;
scatter(distCBDr(empwpl06>0),pccTT_end_HHbar_Cw(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the CBD');
ylabel('Change in wage');
title('Car ban: % change in wage');
legend('Open city','Closed city','Location','best');
%saveas(gcf,'figs/endcft/FIG_TT_OpenClosed_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_TT_OpenClosed_w.png','-dpng','-r300');
% As residents move towards the center and firms move to the suburbs - so
% that the city becomes less specialized - wages fall in the center and
% increase in suburbs due to changes in density. In general the effects are
% more negative in the open city case since the city loses employment and
% agglomeration effects.

% Residential floor space price
pccTT_end_HHbar_CQ=pccTT_end_HHbar(:,5);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(emprsd06>0),pccTT_end_CQ(emprsd06>0));
hold on;
scatter(distCBDr(emprsd06>0),pccTT_end_HHbar_CQ(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the CBD');
ylabel('Change in residential floor space price');
title('Car ban: % residential floor space price');
legend('Open city','Closed city','Location','best');
%saveas(gcf,'figs/endcft/FIG_TT_OpenClosed_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_TT_OpenClosed_Q.png','-dpng','-r300');
% In the close city there is a lot of residential relocation, so prices
% increase in the center. In the open-city case, the city loses
% employment, so there is less demand and the price effect is mostly
% negative.

% Commercial floor space price
pccTT_end_HHbar_Cq=pccTT_end_HHbar(:,5);
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distCBDr(emprsd06>0),pccTT_end_Cq(emprsd06>0));
hold on;
scatter(distCBDr(emprsd06>0),pccTT_end_HHbar_Cq(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the CBD');
ylabel('Change in floor space price');
title('Car ban: % floor space price');
legend('Open city','Closed city','Location','best');
%saveas(gcf,'figs/endcft/FIG_TT_OpenClosed_HR.png');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_TT_OpenClosed_q.png','-dpng','-r300');

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
[Cendog_cft,Cucprob_cft,HH_cft,CUbar_cft,Cconverge]=ussmodendog(param,fund,tt06,nobs06,Ubar06);
wtime=toc;
display('>>>> Time to solve for equilibrium (mins) <<<<');
wtime./60
% Save data for illustration
save('data/output/endogcftual_AE','Cendog_cft','Cucprob_cft','HH_cft','CUbar_cft','Cconverge');

%%% Store results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change in total employment
delta_HH_AE_end = HH_cft - HH                                               % city gains 300k workers 
% Compute percentage changes and map
pccAE_end = (Cendog_cft ./ Cendog -1).*100; 
pccAE_end(pccAE_end>100)=100;   % Truncate changes at -50% and +100% for a clarity of the map
pccAE_end(pccAE_end<-50)=-50;   % Truncate changes at -50% and +100% for a clarity of the map


%%% Inpect results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wages
display('>>>> Mapping change in wage <<<<');
pccAE_end_Cw=pccAE_end(:,1);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_Cw,'East Berlin productivity gain, open city: % change in wage','figs/endcft','MAP_pcc_EastProd_w') 
% Center becomes more expensive

% Residence employment
display('>>>> Mapping change in residence employment <<<<');
pccAE_end_CHR=pccAE_end(:,8);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_CHR,'Car ban, open city: % change in residence employment','figs/endcft','MAP_pcc_EastProd_HR') 
% Remaining residents leave suburbs and move towards the center 

% Workplace employment
display('>>>> Mapping change in workplace employment <<<<');
pccAE_end_CHM=pccAE_end(:,7);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_CHM,'East Berlin productivity gain, open city: % change in workplace employment','figs/endcft','MAP_pcc_EastProd_HM') 
% Firms are croweded out of the center and become relatively more present
% in suburbs - blocks become less specialized

% floorspace prices
display('>>>> Mapping change in floor space prices <<<<');
pccAE_end_CQ=pccAE_end(:,5);
pccAE_end_Cq=pccAE_end(:,6);
RESULT = MAPIT('../shapefile/Berlin4matlab',pccAE_end_CQ,'East Berlin productivity gain, open city: % change in floor space prices','figs/endcft','MAP_pcc_EastProd_Q') 
% Center becomes more expensive


% Save counterfactual results for later comparison
pccAE_end = pccAE_end;
save('data/output/AEcft_end','delta_HH_AE_end','pccAE_end');   

%%% Compare to closed-city case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data/output/AEcft_end_HHbar','UpccAE_end_HHbar','pccAE_end_HHbar'); 

% Workplace employment
pccAE_end_CHM=pccAE_end(:,7);
pccAE_end_CHM_HHbar=pccAE_end_HHbar(:,7);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(empwpl06>0),pccAE_end_CHM(empwpl06>0));
hold on;
scatter(distWALLrunning(empwpl06>0),pccAE_end_CHM_HHbar(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in workplace employment');
title('East Berlin productivitiy upgrade: % change in workplace employment');
legend('Open city','Closed city','Location','best');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_EastProd_EndExog_HM.png','-dpng','-r300');
% Since total employment increases we see generally an upward shift in the
% open city relative to the closed city

% Residence employment
pccAE_end_CHR=pccAE_end(:,9);
pccAE_end_CHR_HHbar=pccAE_end_HHbar(:,9);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(emprsd06>0),pccAE_end_CHR(emprsd06>0));
hold on;
scatter(distWALLrunning(emprsd06>0),pccAE_end_CHR_HHbar(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in residence employment');
title('East Berlin productivitiy upgrade: % change in residence employment');
legend('Open city','Closed city','Location','best');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_EastProd_EndExog_HR.png','-dpng','-r300');
% Residence employment changes are smoothed across the former border due to
% commuting
% Since total employment increases we see generally an upward shift in the
% open city relative to the closed city

% Wages
pccAE_end_Cw=pccAE_end(:,1);
pccAE_end_Cw_HHbar=pccAE_end_HHbar(:,1);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(empwpl06>0),pccAE_end_Cw(empwpl06>0));
hold on;
scatter(distWALLrunning(empwpl06>0),pccAE_end_Cw_HHbar(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in wage');
title('East Berlin productivitiy upgrade: % change in wage');
legend('Open city','Closed city','Location','best');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_EastProd_EndExog_w.png','-dpng','-r300');
% saveas(gcf,'figs/endvsexog/FIG_EastProd_EndExog_w.png');
% With elastic labour supply to the city, the increase in wage that follows
% from a demand shift is smaller. Therefore, the open-city gradient is
% shifted downwards

% Residential floor space prices
pccAE_end_CQ=pccAE_end(:,5);
pccAE_end_CQ_HHbar=pccAE_end_HHbar(:,5);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(emprsd06>0),pccAE_end_CQ(emprsd06>0));
hold on;
scatter(distWALLrunning(emprsd06>0),pccAE_end_CQ_HHbar(emprsd06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in residential floor space price');
title('East Berlin productivitiy upgrade: % change in residential floor space price');
legend('Open city','Closed city','Location','best');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_EastProd_EndExog_Q.png','-dpng','-r300');
% Resembles residence employment, which is a major driver of demand

% Commercial floor space prices
pccAE_end_Cq=pccAE_end(:,6);
pccAE_end_Cq_HHbar=pccAE_end_HHbar(:,6);
distWALLrunning = distWALLr;
index = dummywestr == 1;
distWALLrunning(index) = distWALLr(index) .* -1;
clf;
fig = figure;
set(fig,'Position',[100,100,1200,600]); % [left, bottom, width, height]
scatter(distWALLrunning(empwpl06>0),pccAE_end_Cq(empwpl06>0));
hold on;
scatter(distWALLrunning(empwpl06>0),pccAE_end_Cq_HHbar(empwpl06>0));
line(xlim,[0 0],'Color','k','LineStyle','--');
line([0 0],ylim,'Color','k','LineStyle','--');
hold off;
xlabel('Distance from the Wall');
ylabel('Change in commercial floor space price');
title('East Berlin productivitiy upgrade: % change in commercial floor space price');
legend('Open city','Closed city','Location','best');
set(fig,'Papersize',[6 3]);
print(fig,'figs/openvsclose/FIG_EastProd_EndExog_q.png','-dpng','-r300');
% Quite similar to residential floor space prices, but the decay is
% steeper, steeply decaying agglomeration effects seem to play a role

display('<<< Counterfactuals with endogenous agglomeration and exogenous reservation utility successfully completed >>>>')