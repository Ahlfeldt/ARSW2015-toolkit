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
%%% This version generates a smaller working file that containy only    %%%
%%% 1986 data                                                           %%%
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
%{
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
%}
global alpha beta kappa epsilon lambda delta eta rho;
global fwestd feastd fwestr feastr; 
global ldistCBDdw ldistDIVdw ldistOUTdw ldistCBDrw ldistDIVrw ldistOUTrw;
%global obsvar36 tt36 nobs36 wage_36;
global obsvar86dw nobs86dw wage_86dw; % tt86dw
global obsvar86rw tt86rw nobs86rw wage_86rw;
%global obsvar06 tt06 nobs06 wage_06;

% **********************************;
% **** READ IN CALIBRATION DATA ****;
% **********************************;

%if bigdata==1;
name06='data/caldata2006_ren.csv';
name36='data/caldata1936_div.csv';
name86ren='data/caldata1986_ren.csv';
name86div='data/caldata1986_div.csv';
distname06='data/blockdistmatrix_ren.csv';
distname36='data/blockdistmatrix_div.csv';
ttfin06='data/ttfinal_2006_ren.csv';
ttpub36='data/ttpublic_1936_div.csv';
ttfin86div='data/ttfinal_1986_div.csv';
ttfin86ren='data/ttfinal_1986_ren.csv';
%end;
% Read in 06 Berlin data;
%data06 = csvread(name06);
% Read in 86ren Berlin data;
data86r = csvread(name86ren);
% Travel time 06 matrix;
%tt06 = csvread(ttfin06);
%tt06(:,1)=[];
% Travel time 86 reunification matrix;
tt86rw = csvread(ttfin86ren);
tt86rw(:,1)=[];
% Distance matrix 06;
%distvec06 = csvread(distname06);
%distvec06 = distvec06(:,3);
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
% **** PREPARE 1936 DATA ****;
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
%{
nobs36=size(data36,1);
block36=data36(:,1);
areaidd=data36(:,2);
bzk36=data36(:,3);
dummywestd=data36(:,5);
latituded=data36(:,6);
longituded=data36(:,7);
floor36=data36(:,8);
gfz36=data36(:,9);
area36=data36(:,10);
empwpl36=data36(:,12);
emprsd36=data36(:,13);
distCBDd=data36(:,14);
distWALLd=data36(:,15);
distOUTd=data36(:,16);
USB36d=data36(:,17);
USB86d=data36(:,18);
distTd=data36(:,19);
distKUd=data36(:,20);
Xd=data36(:,21);
Yd=data36(:,22);

% Create distance matrix from vector;
[dist36] = reshape(distvec36,nobs36,nobs36);
% Create measure of internal distance 36;
intdist36=(2/3).*((area36/pi).^(1/2));
x=1;
while x<=nobs36;
     dist36(x,x)=intdist36(x);
     x=x+1;
end;

% Scale employment workplace so total employment workplace equals;
% total employment residence;
[empwpl36]=empadjust(empwpl36,emprsd36);
totemp36=sum(empwpl36);
test=sum(emprsd36);
%}
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
%{
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
%}
% *******************************;
% **** WEST BERLIN VARIABLES ****;
% *******************************;
%{
% Delete East Berlin from the 36 data;
fwestd=find(dummywestd==1);
feastd=find(dummywestd==0);
[areaiddw]=eastdel(areaidd,feastd);
[floor36w]=eastdel(floor36,feastd);
[area36w]=eastdel(area36,feastd);
[empwpl36w]=eastdel(empwpl36,feastd);
[emprsd36w]=eastdel(emprsd36,feastd);
area86dw=area36w;
[distCBDdw]=eastdel(distCBDd,feastd);
[distOUTdw]=eastdel(distOUTd,feastd);
[distWALLdw]=eastdel(distWALLd,feastd);
[distTdw]=eastdel(distTd,feastd);
[distKUdw]=eastdel(distKUd,feastd);
[dummywestdw]=eastdel(dummywestd,feastd);
[latitudedw]=eastdel(latituded,feastd);
[longitudedw]=eastdel(longituded,feastd);
[gfz36dw]=eastdel(gfz36,feastd);
[USB36dw]=eastdel(USB36d,feastd);
[USB86dw]=eastdel(USB86d,feastd);
[Xdw]=eastdel(Xd,feastd);
[Ydw]=eastdel(Yd,feastd);

% Delete East Berlin from the 06 data;
fwestr=find(dummywestr==1);
feastr=find(dummywestr==0);
[areaidrw]=eastdel(areaidr,feastr);
[floor06w]=eastdel(floor06,feastr);
[area06w]=eastdel(area06,feastr);
[empwpl06w]=eastdel(empwpl06,feastr);
[emprsd06w]=eastdel(emprsd06,feastr);
[area86rw]=area06w;
[distCBDrw]=eastdel(distCBDr,feastr);
[distOUTrw]=eastdel(distOUTr,feastr);
[distWALLrw]=eastdel(distWALLr,feastr);
[distTrw]=eastdel(distTr,feastr);
[distKUrw]=eastdel(distKUr,feastr);
[dummywestrw]=eastdel(dummywestr,feastr);
[latituderw]=eastdel(latituder,feastr);
[longituderw]=eastdel(longituder,feastr);
[gfz06rw]=eastdel(gfz06,feastr);
[USB36rw]=eastdel(USB36r,feastr);
[USB86rw]=eastdel(USB86r,feastr);
[Xrw]=eastdel(Xr,feastr);
[Yrw]=eastdel(Yr,feastr);

% Extract 1986 bilateral data;
[dist86dw]=eastddel(dist36,feastd);
% Extract 1986 bilateral data;
[dist86rw]=eastddel(dist06,feastr);
%}
% Extract 1986 data;
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
data06 = csvread(name06);
dummywestr=data06(:,5);
clear data06
fwestr=find(dummywestr==1);
%[data86dw]=eastdel(data86d,feastd);
[data86rw]=eastdel(data86r,feastr);
%nobs86dw=size(data86dw,1);
nobs86rw=size(data86rw,1);
%block86dw=data86dw(:,1);
block86rw=data86rw(:,1);
%bzk86dw=data86dw(:,3);
bzk86rw=data86rw(:,3);
%floor86dw=data86dw(:,8);
floor86rw=data86rw(:,8);
%gfz86dw=data86dw(:,9);
gfz86rw=data86rw(:,9);
%empwpl86dw=data86dw(:,12);
empwpl86rw=data86rw(:,12);
%emprsd86dw=data86dw(:,13);
emprsd86rw=data86rw(:,13);

% Scale employment workplace so total employment workplace equals;
% total employment residence;
% 1936 West;
%[empwpl36w]=empadjust(empwpl36w,emprsd36w);
%totemp36w=sum(empwpl36w);
% 2006 West;
%[empwpl06w]=empadjust(empwpl06w,emprsd06w);
%totemp06w=sum(empwpl06w);
% 1986 West division;
%[empwpl86dw]=empadjust(empwpl86dw,emprsd86dw);
%totemp86dw=sum(empwpl86dw);
% 1986 West reunification;
[empwpl86rw]=empadjust(empwpl86rw,emprsd86rw);
totemp86rw=sum(empwpl86rw);

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

%if bigdata==1;
save('data/input/prepdata_big_TD86');
%end;

display('>>>> File Completed Successfully <<<<');
