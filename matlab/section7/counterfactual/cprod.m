%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program is part of the replication directory                   %%%
%%% This version has been commented                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This function decomposes prodocutivities A into an endogenous component
%%% that depends on density and a fundamental exogenous component a

function [a,Ups] = cprod(obsdata,distvar,noj,A);                            % wage has been removed as an argument since it is redundant
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n    
        % A are adjusted producitivities recovered by cmodexorg.m
    % This program produces the following outputs
        % a is the exogenous component in adjusted productivities
        % Ups is the endogenous component
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen

% Declaring scalars as globals so they can be accessed and modified by any function    
global alpha beta epsilon kappa lambda delta rho eta;

% Extracting variables from the four key input variables
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area
  
% Only take dimensions with positive employment
Iwpl = (HMT~=0); Inwpl = (HMT==0);                                          % Two variables that define when workplace employment is positive or not positive (zero)
Irsd = (HRT~=0); Inrsd = (HRT==0);                                          % Two variables that define when residence employment is positive or not positive (zero)
nto = sum(Iwpl);                                                            % The number of observations with positive workplace employment
nfrom = sum(Irsd);                                                          % The number of observations with positive residence employment
EHMT = HMT(Iwpl);                                                           % Generates a vector of workplace employment only containing observations with positive workplace emplyoment (it has nto observations)
EHRT = HRT(Irsd);                                                           % Generates a vector of residence employment only containing observations with positive residence emplyoment (it has nfrom observations)
EKM = K(Iwpl);                                                              % Generates a vector of geographic area only only containing observations with positive workplace employment
% Ewage = wage(Iwpl);
EA = A(Iwpl);                                                               % Generates a vector of adjusted productivities only containing obserations with positive workplace employment
HH=sum(EHMT);                                                               % Compute total workplace employment

% We compute UPSILON using equation (20);
dd=exp(-delta.*distvar(Iwpl,Iwpl));                                         % We create a nto x nto matix of spatial weights, the exp^(-delta*tau_js) component in Equation (20); notice that that rows and columns that are not selected by the index (Iwpl,Iwpl) are residential and have zero workplace employment. Thefore, these can be ignored.
EUps=dd*(EHMT./EKM);                                                        % We compute an nto x 1 vector of employment densities, the j-specific component in the summation in Equation (20)
Ups=zeros(noj,1);                                                           % We compute UPSILON by way of matrix multiplication of the nto x nto matrix of spatial weights with the nto x 1 vector of employment densities. The matrix multiplication sums employment over all j weighted by the spatial weight.    
Ups(Iwpl)=EUps;                                                             % Create N x 1 vector of Upsilon. Places without workplace employment obtain theory-consistent zeros.

% PRODUCTIVITY;
a=zeros(noj,1);                                                             % Placeholder for fundamental productivities
a(Iwpl)=EA./(EUps.^lambda);                                                 % Create N x 1 vector of fundamental productities. Places without workplace employment obtain theory-consistent zeros.


