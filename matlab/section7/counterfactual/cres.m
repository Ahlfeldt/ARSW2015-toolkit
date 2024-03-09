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

%%% This function decomposes amenities B into an endogenous component
%%% that depends on density and a fundamental exogenous component b

function [b,Ome] = cres(obsdata,distvar,noj,B)                              % wage has been removed as an argument since it is redundant
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n    
        % B are adjusted amenities recovered by cmodexorg.m
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
%Ewage = wage(Iwpl);
EKM = K(Iwpl);                                                              % Generates a vector of geographic area only only containing observations with positive workplace employment
EKR = K(Irsd);                                                              % Generates a vector of geographic area only only containing observations with positive residence employment
EB = B(Irsd);                                                               % Generates a vector of adjusted amenities only containing obserations with positive residence employment
HH=sum(EHMT);                                                               % Compute total workplace employment

% We compute OMEGA using equation (21);
rr=exp(-rho.*distvar(Irsd,Irsd));                                           % We create a nfrom x nfrom matix of spatial weights, the exp^(-rho*tau_js) component in Equation (21); notice that that rows and columns that are not selected by the index (Irsd,Irsd) are commercial and have zero residence employment. Thefore, these can be ignored.
EOme=rr*(EHRT./EKR);                                                        % We compute an nfrom x 1 vector of residence employment densities, the j-specific component in the summation in Equation (21)
Ome=zeros(noj,1);                                                           % We compute OMEGA by way of matrix multiplication of the nfrom x nfrom matrix of spatial weights with the nfrom x 1 vector of employment densities. The matrix multiplication sums employment over all j weighted by the spatial weight.    
Ome(Irsd)=EOme;                                                             % Create N x 1 vector of Omega. Places without residence employment obtain theory-consistent zeros.

% RESIDENTIAL FUNDAMENTALS;
b=zeros(noj,1) ;                                                            % Placeholder for fundamental amenities
b(Irsd)=EB./(EOme.^eta);                                                    % Create N x 1 vector of fundamental amenities. Places without residence employment obtain theory-consistent zeros.



