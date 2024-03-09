%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The orginal version of this file was part of the replication        %%%
%%% directory to Ahlfeldt, Redding, Sturm, Wolf (2015)                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This version has modified by Gabriel M. Ahlfeldt in 2024            %%%
%%% It produces a smaller working file that containy only 2006 data     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
clf;
colormap default; 
format bank;
close all;

% User==1 is Steve Windows desktop;
% User==2 is Steve Mac Pro desktop;
% User==5 is Berlin Humboldt server;

% bigdata==1 is the full dataset;

user=1;

% Select data;
bigdata=1;

if user==1;
    cd 'E:/Data/Gabriel/Ahfleldt/_QSE/ARSW/FinalBerlin/matlab';
elseif user==2;
    cd '/Users/reddings/Dropbox/SRDS/FinalBerlin/matlab/';
elseif user==5;
    cd '/work/eodwiwi/FinalBerlin/matlab/';
end;

global alpha beta kappa epsilon lambda delta eta rho;
global fwestd feastd fwestr feastr; 
global ldistCBDdw ldistDIVdw ldistOUTdw ldistCBDrw ldistDIVrw ldistOUTrw;
global obsvar06 tt06 nobs06 wage_06;

% **********************************;
% **** READ IN CALIBRATION DATA ****;
% **********************************;

if bigdata==1;
name06='data/caldata2006_ren.csv';
%name36='data/caldata1936_div.csv';
%name86ren='data/caldata1986_ren.csv';
%name86div='data/caldata1986_div.csv';
distname06='data/blockdistmatrix_ren.csv';
%distname36='data/blockdistmatrix_div.csv';
ttfin06='data/ttfinal_2006_ren.csv';
%ttpub36='data/ttpublic_1936_div.csv';
%ttfin86div='data/ttfinal_1986_div.csv';
%ttfin86ren='data/ttfinal_1986_ren.csv';
end;
% Read in 06 Berlin data;
data06 = csvread(name06);
% Read in 86ren Berlin data;
%data86r = csvread(name86ren);
% Travel time 06 matrix;
tt06 = csvread(ttfin06);
tt06(:,1)=[];
% Travel time 86 reunification matrix;
%tt86rw = csvread(ttfin86ren);
%tt86rw(:,1)=[];
% Distance matrix 06;
distvec06 = csvread(distname06);
distvec06 = distvec06(:,3);
% Read in 36 Berlin data;
%data36 = csvread(name36);
% Read in 86div Berlin data;
%data86d = csvread(name86div);
% Travel time 36 matrix;
%tt36 = csvread(ttpub36);
%tt36(:,1)=[];
% Travel time 86 division matrix;
%tt86dw = csvread(ttfin86div);
%tt86dw(:,1)=[];
% Distance matrix 36;
%distvec36 = csvread(distname36);
%distvec36 = distvec36(:,3);


% ***************************;
% **** PREPARE 2006 DATA ****;
% ***************************;

% Column 1 : block_id
% Column 2 : area_id 
% Column 3 : bezirk1937 
% Column 4 : year 
% Column 5 : dummywest
% Column 6 : latitude 
% Column 7 : longitude 
% Column 8 : floor 
% Column 9 : gfz 
% Column 10 : finhalt 
% Column 11 : factories
% Column 12 : emp_wpl 
% Column 13 : emp_rsd 
% Column 14 : distCBD 
% Column 15 : d_wall 
% Column 16 : d_outb
% Column 17 : USB36;
% Column 18 : USB86;
% Column 19 : dt1936;
% Column 20 : distku;
% Column 21 : x_coord;
% Column 22 : y_coord;

nobs06=size(data06,1);
block06=data06(:,1);
areaidr=data06(:,2);
bzk06=data06(:,3);
dummywestr=data06(:,5);
latituder=data06(:,6);
longituder=data06(:,7);
floor06=data06(:,8);
gfz06=data06(:,9);
area06=data06(:,10);
empwpl06=data06(:,12);
emprsd06=data06(:,13);
distCBDr=data06(:,14);
distWALLr=data06(:,15);
distOUTr=data06(:,16);
USB36r=data06(:,17);
USB86r=data06(:,18);
distTr=data06(:,19);
distKUr=data06(:,20);
Xr=data06(:,21);
Yr=data06(:,22);

% Create distance matrix from vector;
[dist06] = reshape(distvec06,nobs06,nobs06);
% Create measure of internal distance;
intdist06=(2/3).*((area06/pi).^(1/2));
x=1;
while x<=nobs06;
     dist06(x,x)=intdist06(x);
     x=x+1;
end;

% Scale employment workplace so total employment workplace equals;
% total employment residence;
[empwpl06]=empadjust(empwpl06,emprsd06);
totemp06=sum(empwpl06);
test=sum(emprsd06);

% *****************************;
% **** X and Y coordinates ****;
% *****************************;
%{
% Matrices of differences in Cartesian coordinates;
% Division;
Xd1=repmat(Xd,1,nobs36);
Xd2=repmat(Xd',nobs36,1);
Yd1=repmat(Yd,1,nobs36);
Yd2=repmat(Yd',nobs36,1);
Xdmat=Xd1-Xd2;
Ydmat=Yd1-Yd2;
% Reunification
Xr1=repmat(Xr,1,nobs06);
Xr2=repmat(Xr',nobs06,1);
Yr1=repmat(Yr,1,nobs06);
Yr2=repmat(Yr',nobs06,1);
Xrmat=Xr1-Xr2;
Yrmat=Yr1-Yr2;
% Division West;
Xdw1=repmat(Xdw,1,nobs86dw);
Xdw2=repmat(Xdw',nobs86dw,1);
Ydw1=repmat(Ydw,1,nobs86dw);
Ydw2=repmat(Ydw',nobs86dw,1);
Xdwmat=Xdw1-Xdw2;
Ydwmat=Ydw1-Ydw2;
% Reunification West;
Xrw1=repmat(Xrw,1,nobs86rw);
Xrw2=repmat(Xrw',nobs86rw,1);
Yrw1=repmat(Yrw,1,nobs86rw);
Yrw2=repmat(Yrw',nobs86rw,1);
Xrwmat=Xrw1-Xrw2;
Yrwmat=Yrw1-Yrw2;
% Absolute differences;
Xdwmat=abs(Xdwmat);
Ydwmat=abs(Ydwmat);
Xrwmat=abs(Xrwmat);
Yrwmat=abs(Yrwmat);
% Map to kilometers;
Xdwmat=Xdwmat./1000;
Ydwmat=Ydwmat./1000;
Xrwmat=Xrwmat./1000;
Yrwmat=Yrwmat./1000;
% Clear variables;
clear longituded latituded longituder latituder;
clear coslatituded coslongituded sinlatituded sinlongituded;
clear coslatituder coslongituder sinlatituder sinlongituder;
clear a_ f_ e_2 tempd tempd normlatd tempr tempr normlatr;
clear Xd1 Xd2 Yd1 Yd2 Xr1 Xr2 Yr1 Yr2;
clear Xdw1 Xdw2 Ydw1 Ydw2 Xrw1 Xrw2 Yrw1 Yrw2;
%}
% *******************;
% **** Save Data ****;
% *******************;

clear dist06 distvec06

if bigdata==1;
save('data/input/prepdata_big_TD');
end;

display('>>>> File Completed Successfully <<<<');
