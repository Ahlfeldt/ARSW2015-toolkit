%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by SJR, 02/2015                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This version of this program is part of the replication directory   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outvar] = modbezirk(invar,noj)

% ************************;
% **** Modern Bezirke ****;
% ************************;

% Old Bezirke;
% 1 :  Charlottenburg;
% 2 :  Friedrichshain;
% 3 : Hellersdorf;
% 4 : Hohenschönhausen;
% 5 : Kreuzberg;
% 6 : Köpenick;
% 7 : Lichtenberg;
% 8 : Marzahn;
% 9 : Mitte;
% 10 : Neukölln;
% 11 : Pankow;
% 12 : Prenzlauer Berg;
% 13 : Reinickendorf;
% 14 : Schöneberg;
% 15 : Spandau;
% 16 : Steglitz;
% 17 : Tempelhof;
% 18 : Tiergarten;
% 19 : Treptow;
% 20 : Wedding;
% 21 : Weißensee;
% 22 : Wilmersdorf;
% 23 : Zehlendorf;

% Modern bezirke;
% 1 : Charlottenburg-Wilmersdorf;
% 2 : Friedrichshain-Kreuzberg;
% 3 : Lichtenberg;
% 4 : Marzahn-Hellersdorf;
% 5 : Mitte;
% 6 : Neukolln;
% 7 : Pankow;
% 8 : Reinickendorf;
% 9 : Spandau;
% 10 : Steglitz-Zehlendorf;
% 11 : Tempelhof-Schoneberg;
% 12 : Treptow-Kopenick;
% Modern bezirke;
outvar=zeros(noj,1);
% 1 : Charlottenburg-Wilmersdorf;
outvar(invar==1)=1;
outvar(invar==22)=1;
% 2 : Friedrichshain-Kreuzberg;
outvar(invar==2)=2;
outvar(invar==5)=2;
% 3 : Lichtenberg;
outvar(invar==7)=3;
outvar(invar==4)=3;
% 4 : Marzahn-Hellersdorf;
outvar(invar==8)=4;
outvar(invar==3)=4;
% 5 : Mitte;
outvar(invar==9)=5;
outvar(invar==18)=5;
outvar(invar==20)=5;
% 6 : Neukolln;
outvar(invar==10)=6;
% 7 : Pankow;
outvar(invar==11)=7;
outvar(invar==12)=7;
outvar(invar==21)=7;
% 8 : Reinickendorf;
outvar(invar==13)=8;
% 9 : Spandau;
outvar(invar==15)=9;
% 10 : Steglitz-Zehlendorf;
outvar(invar==16)=10;
outvar(invar==23)=10;
% 11 : Tempelhof-Schoneberg;
outvar(invar==17)=11;
outvar(invar==14)=11;
% 12 : Treptow-Kopenick;
outvar(invar==19)=12;
outvar(invar==6)=12;
