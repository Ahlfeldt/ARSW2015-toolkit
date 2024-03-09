%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: GMA, 03/2024                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program is not part of the replication directory               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This funciton computes the reservation utility level in the city
function [Ubar] = ubar(B,param,fund,distvar,noj);
       % This program uses the following inputs
       % B06 is the adjusted amenity vector recovered by cmodexog.m
       % param is a 1 x k vector of parameters epsi and kappaepsi
        % fund is a 1 x k vector of N x 1 vectors of exogenous fundamentals
            % and starting values of endogenous variables
                % Productivities, amenities, density of development, area,
                % floor space prices, workplace employment, residence
                % employment, commercial floor space, residential floor
                % space, total floor space, adjusted wages, total worker
                % income      
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
    % This program produces the following outputs
        % Ubar is the exogenous tility level.
    % The names of the inputs need to correspond to objects that exist in
        % the workspace. The names of the outputs can be freely chosen

% Declaring scalars as globals so they can be accessed and modified by any function
global alpha beta  

% Read parameters from parameter input vector
epsilon = param(2);
kappa = param(1);

% Read fundamentals and observed variables from variable input vector
a = fund(:,1);                                                              % Fundamental productivitiy
b = fund(:,2);                                                              % Fundamental amenity
floor = fund(:,5);                                                          % Floor space prices
wage = fund(:,8);     
% Wages
% Generate index variables
Ito=(a~=0);                                                                 % Index all places with positive fundamental productivitiy                                                           
Ifrom=(b~=0);                                                               % Index all places with positive fundamental amenity
nto=sum(Ito);                                                               % Number of places with positive fundamental productivity
nfrom=sum(Ifrom);                                                           % Number of places with positive fundamental amenity

% Compute components of Eq. (9)
gammaf=gamma((epsilon-1)./epsilon);                                         % Gamma from Eq. (9) defining expected utility and Ubar
d_ij_eps=exp(-epsilon.*kappa.*distvar(Ifrom,Ito));                          % Commuting cost component, notice that d_ij = exp(-kappa*tau_ij)
EQQ=(floor(Ifrom).^(-(1-beta).*epsilon))*ones(1,nto);                       % Compute floor space price component in Eq. (9)
EBB=(B(Ifrom)*ones(1,nto)).^epsilon;                                        % Compute amenity component in Eq. (9)
EWW=((wage(Ito)*ones(1,nfrom))').^epsilon;                                  % Compute wage component in Eq. (9)
Ephi_ij=d_ij_eps.*EBB.*EQQ.*EWW;                                            % Multiply the above to compute Phi_ij defined in Eq. (4)
Ephi=sum(sum(Ephi_ij));                                                     % Compute Phi from Eq. (4)
Ubar=gammaf.*(Ephi.^(1./epsilon))                                           % Compute reservation utility level according to Eq. (9)

display('<<<<<<<<<<<<< Reservation utility U_bar computed >>>>>>>>>>>>>>>')
