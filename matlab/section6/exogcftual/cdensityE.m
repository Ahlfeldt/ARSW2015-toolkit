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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vout,LD,LM,LR,thetaout] = cdensityE(obsdata,A,wage,vv,noj);
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % A is adjusted productivities, solved by comegaoptC.m
        % Wage is the vector of adjusted wages solved by either
            % cmodexog.m
        % vv is total expected worker income solved by  cmodexog.m
        % noj is a scalar that defines the number of observations n     
    % This program produces the following outputs
        % Vout is the density of development varphi
        % LD is total requilibrium floor space
        % LM id commercial equilibrium floor space
        % LR is residential equilibrium floor space
        % thetaout is the share of residential floor space at total floor
            % space
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen
    
% Declaring scalars as globals so they can be accessed and modified by any function    
global alpha beta;

% Extracting variables from the four key input variables
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area

% Only take dimensions with positive employment
Iwpl = (HMT~=0); Inwpl = (HMT==0);                                          % Two variables that define when workplace employment is positive or not positive (zero)
Irsd = (HRT~=0); Inrsd = (HRT==0);                                          % Two variables that define when residence employment is positive or not positive (zero)

% *********************;
% **** Floor space ****;
% *********************;

% Commercial land demand;
LM=zeros(noj,1);                                                            % We generate a vector for commercial floor space demand (1-theta)L
LM(Iwpl)=((((1-alpha).*A(Iwpl))./QT(Iwpl)).^(1./alpha)).*HMT(Iwpl);         % We generate commercial floor space (1-theta)L using equation S.30 in the supplement. We index to change values for locations with positive workplace employment. The other locations remain theory-consistent zero.
                                                            
% Residential floor space demand;
LR=zeros(noj,1);                                                            % We generate a vector for residential floor space demand (theta)L
LR(Irsd)=((1-beta).*vv(Irsd))./QT(Irsd);                                    % We generate residential floor space (theta)L using equation S.29 in the supplement. We index to change values for locations with positive residence employment. The other locations remain theory-consistent zero.

% Total floor space demand;
LD=LM+LR;                                                                   % Total floor space demand is simply the sum of commercial and residential floor space demand                                                              
Vout=LD./(K.^0.75);                                                         % We solve for density of development varphi using Eq. S31
thetaout=LM./LD;                                                            % Share of commercial floor space theta is residential floor space
