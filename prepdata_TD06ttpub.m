%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The orginal version of this file was part of the replication        %%%
%%% directory to Ahlfeldt, Redding, Sturm, Wolf (2015)                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This converts the 2006 public transport matrix into a matlab file   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
clf;
colormap default; 
format bank;
close all;


ttpub06='data/ttpublic_2006_ren.csv';
ttpub06 = csvread(ttpub06);
ttpub06(:,1)=[];

% *******************;
% **** Save Data ****;
% *******************;

save('data/ttpublic_2006_ren');


display('>>>> File Completed Successfully <<<<');
