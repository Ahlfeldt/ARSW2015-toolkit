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

%%% This function adjust the levels adjusted produtivities and adjusted
%%% amenities recovered from the sequential procedure called by calcal_TD
%%% so that when solving the model obtain the correct total
%%% population. Using this program is important when conducting
%%% counterfactual analyses starting from the productivities and amenities
%%% recovered by the sequential procedure. 
%%%
%%% Intuitively, we rescale wages so that they are consistent with
%%% productivities that have a geometric mean of one and, conditional on 
%%% rescaled wages, we amenities so that we match the city population in
%%% the model (Phi = city population)

function [Aout,Bout,wageout] = calcal_adj_TD(obsdata,distvar,noj,A,B)
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
        % A contains guesses of producitivities
        % B contains guesses of amenities
    % This program produces the following outputs
        % Aout is is a vector of updated inverted adjusted productivities
        % Bout is is a vector of updated inverted adjusted amenities
    
% Declaring scalars as globals so they can be accessed and modified by any function    
    global alpha beta kappa epsilon;        
        
% Read key variables from obsdata input
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area

% Prepare variables for blocks with postive reisdence or workplace amenities    
    Iwpl = (HMT~=0); Inwpl = (HMT==0);                                      % Two variables that define when workplace employment is positive or not positive (zero)
    Irsd = (HRT~=0); Inrsd = (HRT==0);                                      % Two variables that define when residence employment is positive or not positive (zero)
    nwpl = sum(Iwpl);                                                       % The number of observations with positive workplace employment
    nrsd = sum(Irsd);                                                       % The number of observations with positive residence employment 
    EHMT = HMT(Iwpl);                                                       % Generates a vector of workplace employment only containing observations with positive workplace emplyoment (it has nto observations)
    EHRT = HRT(Irsd);                                                       % Generates a vector of residence employment only containing observations with positive residence emplyoment (it has nfrom observations)
    wage = zeros(noj,1);                                                    % Placeholder for wage vector to be generated 
    Ewage = wage(Iwpl);                                                     % Generates a vector of guesses of wages only containing observations with positive workplace employment
    EA = A(Iwpl);                                                         % Generates a vector of adjusted productivities only containing obserations with positive workplace employment
    EA=EA./geomean(EA);                                                     % Normalize productivities to have a geometric mean of one
                                                                            % This is important since in the sequential procedure we have normalized adjusted wages to geometric mean of one which does not imply that productivities have a geometric mean of one        
    EB = B(Irsd);                                                         % Generates a vector of adjusted amenities only containing obserations with positive residence employment
    EKM = K(Iwpl);                                                          % Generates a vector of geographic area only only containing observations with positive workplace employment
    EKR = K(Irsd);                                                          % Generates a vector of geographic area only only containing observations with positive residence employment    
    HH=sum(EHMT);   

% Update Ewage to reflect normalized productivities
Ewage=(((1-alpha)./QT(Iwpl)).^((1-alpha)./alpha)).*alpha.*(EA.^(1./alpha)); % We use Equation (12) which combines first-order condition and zero-profit conditions.

% Prepare distance variables
Edistvar=distvar(Irsd,Iwpl);                                                % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom osbervations)
clear distvar                                                               % Notice that this time, we are having residence in rows and workplace in columsn, so we have an nfrom x nto matrix
d_ij_eps=exp(-epsilon.*kappa.*Edistvar);                                    % Generates a spatial weight matrix with weights that exponentially falls in bilateral travel time between i and j at rate epsilon*kappa
clear Edistvar   

% Generate population in the model (Phi)
EQQ=repmat(QT(Irsd),1,nwpl);                                                % We assign residence floor space prices to all bilaterals. We take the vector of floor space prices, generate a new vector only for observations with positive residence employment, and replicates it nwpl times to generate a nrsd x nwpl matrix
EQQ=EQQ.^(-(1-beta).*epsilon);                                              % Generate the first component of bilateral commuting probabilities in Eq. 4: floor space prices to the power of (1-beta) x -epsilon 
EBB=repmat(EB,1,nwpl);                                                      % We assign residence amenities to all bilaterals. Takes the vector of amenities for observations with positive residence emplyoment and replicates it to generate an nrds x nwpl matrix
EBB=EBB.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities: B to the power of epsilon 
EWW=repmat(Ewage',nrsd,1);                                                  % We assign workplace wages to all bilaterals. We take the nwpl x 1 vector of adjusted wages, transpose it into an 1 x nwpl vector and replicate in nrsd rows to obtain an nrsd x nwpl matrix
EWW=EWW.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities, wage to the power of epsilon
Ephi_ij=d_ij_eps.*EBB.*EQQ.*EWW;                                            % Now we use the above components to compute the numerator of bilateral commuting probabilities in Eq. 4, phi_ij
clear EQQ EWW EBB                                                           % Clear objects to save memory
Ephi=sum(sum(Ephi_ij));                                                     % We sum over all bilaterals to obtain the denominator in Eq. 4. sum(Ephi_ij) retuns column sums, so we sum over the column sum to get the sum over all bilaterals

% Adjust amenities to match the right population
EB=((HH./Ephi).^(1./epsilon)).*EB;                                          % If the ratio of the total population in the data over the population in the model is positive, we inflate amenities
                                                                            % From normalization of U_bar explained on page 18 in the supplement we have that H = Phi (the term in Brackets in the H equation). This equation reveals that we can move a multiplicative component Badjustment^epsilon in B_i out of the summation and use it to reach any population level H. 
                                                                            % If we have H = "adjustment factor"^epsilon * Phi, the it follows that "adjustment factor" = (H/Phi)^(1/epsilon)
% Compute recaled productivities and amenities for all blocks
Aout=zeros(noj,1);
Aout(Iwpl)=EA;                                                              % Fill locations with positive workplace employment. Those without receive theory-consistent zeros
Bout=zeros(noj,1); 
Bout(Irsd)=EB;                                                              % Fill locations with positive residence employment. Those without receive theory-consistent zeros
wageout=zeros(noj,1); 
wageout(Iwpl)=Ewage ;                                                          % Output updated wage and fill in zero values for locations without workplace employment
display('>>>> Productivities and amenities updated <<<<');

        
       