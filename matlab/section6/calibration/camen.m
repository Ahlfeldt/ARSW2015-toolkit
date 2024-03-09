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

%%% This function uses residence employment, floor space prices and
%%% commuting market access to solve for adjusted amenities using Eq. S.47

function [B,CMA,HRS] = camen(obsdata,distvar,noj,wage,cprob)
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
        % wage is vector of guesses of adjusted wages
        % cprob is the matrix of condional commuting probabilities (rows
            % are workplaces, columns are residences)
    % This program produces the following outputs
        % B is the inverted vector of adjusted amenities 
        % CMA is (residential) commuting market access
        % HRS is share of blocks at residence employment
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen
    
% Declaring scalars as globals so they can be accessed and modified by any function
global alpha beta kappa epsilon;

% Extracting variables from the four key input variables
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area

% Only take dimensions with positive employment;
Iwpl = (HMT~=0); Inwpl = (HMT==0);                                          % Two variables that define when workplace employment is positive or not positive (zero)
Irsd = (HRT~=0); Inrsd = (HRT==0);                                          % Two variables that define when residence employment is positive or not positive (zero)
nto = sum(Iwpl);                                                            % The number of observations with positive workplace employment
nfrom = sum(Irsd);                                                          % The number of observations with positive residence employment 
EHMT = HMT(Iwpl);                                                           % Generates a vector of workplace employment only containing observations with positive workplace emplyoment (it has nto observations)
EHRT = HRT(Irsd);                                                           % Generates a vector of residence employment only containing observations with positive residence emplyoment (it has nfrom observations)
Ewage = wage(Iwpl);                                                         % Generates a vector of adjusted wages only only containing observations with positive workplace employment
EKM = K(Iwpl);                                                              % Generates a vector of geographic area only only containing observations with positive workplace employment                                                       
EKR = K(Irsd);                                                              % Generates a vector of geographic area only only containing observations with positive residence employment    
HH=sum(EHMT);                                                               % Compute total workplace employment
Edistvar=distvar(Irsd,Iwpl);                                                % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom osbervations)
                                                                            % Notice that this time, we are having residence in rows and workplace in columsn, so we have an nfrom x nto matrix
d_ij_eps=exp(-epsilon.*kappa.*Edistvar);                                    % Generates a spatial weight matrix with weights that exponentially falls in bilateral travel time between i and j at rate epsilon*kappa
clear distvar                                                               % clear variable to save memory 

% RESIDENTS SHARE RELATIVE TO GEOMETRIC MEAN;
EHRS=EHRT./sum(EHRT);                                                       % Share of blocks at total residence employment
EHRS=EHRS./geomean(EHRS);                                                   % Normalize by geometric mean

% FLOOR PRICES RELATIVE TO GEOMETRIC MEAN;
EQT=QT(Irsd);                                                               % Floor space prices for the subset of blocks with positive residence employment
EQT=EQT./geomean(EQT);                                                      % Normalize by geometric mean

% COMMUTING ACCESS RELATIVE TO GEOMETRIC MEAN;
ECMA=d_ij_eps*(Ewage.^epsilon);                                             % We perform a matrix multiplication of the spatial weights matrix of dimension nfrom x nto and the nto x 1 vector of adjusted wages. This gives us the nfrom x 1 vector of (residential) commuting market access defined in Section S.3.1.2. 
ECMA=ECMA./geomean(ECMA);                                                   % Normalize commuter market acces by geometric mean

% AMENITIES RELATIVE TO GEOMETRIC MEAN;
EB=(EHRS.^(1./epsilon)).*(EQT.^(1-beta))./(ECMA.^(1./epsilon));             % Now we are ready to invert adjusted amenities using equation S.47
B=zeros(noj,1);                                                             % Create a column vector of length n
B(Irsd)=EB;                                                                 % Map adjusted amenities from the vector containing locations with postive residence employment into the vector containing all observations. Observations with zero residence employment receive a theory-consistent zero amentiy value. 

% COMMUTING ACCESS;
ECMA=d_ij_eps*(Ewage.^epsilon);                                             % We recreate residential commuter market access since we have normalized it before
CMA=zeros(noj,1);                                                           % We create a column vector of length n
CMA(Irsd)=ECMA;                                                             % Map commuter market access from the vector containing locations with postive residence employment into the vector containing all observations. Observations with zero residence employment receive a theory-consistent zero value. 

% RESIDENTIAL EMPLOYMENT SHARE;
EHRS=EHRT./sum(EHRT);                                                       % Finally, compute the share of blocks at total residence employment and map it into a vector containing all locations (those zero residence employment naturally receive a zero share)
HRS=zeros(noj,1);                                                            
HRS(Irsd)=EHRS;




