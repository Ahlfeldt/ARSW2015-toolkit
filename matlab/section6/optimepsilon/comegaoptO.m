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
%%% This version has been commented                                     %%%
%%% It has been modified to be more memory efficient                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This function solves for transformed wages (omega) for given values of 
%%% workplace employment, residence employment, and bilateral commuting 
%%% costs

function [omout,cprob,wconverge,HMC,gap] = comegaoptO(obsdata,distvar,noj,omega);
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, workplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
        % omega is a vector of guesses of transformed wages omega
    % This program produces the following outputs
        % omout is the solved equilibrium vector of trasnformed wages 
        % cprob is the matrix of conditional commuting probabilities (rows
            % are workplaces, columns are residences)
        % wconverge is a dummy indicating that the solver has converged
        % HMC is predicted workplace employment (should be the same as
            % observed workplace employment in case of perfect convergence)
        % gap is the largest difference between any pair of 
            % workplace employment and predicted workplace employment
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen.
        
    % Notice that this function does not really use floor space and area in
    % the obsdata input

% Declaring scalars as globals so they can be accessed and modified by any function
global alpha beta kappaeps;

% Setting the convergence indicator to zero (if convergence is reached it
% will be replaced with 1
wconverge=0;

% Extracting variables from the four key input variables
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area

% Only take dimensions with positive employment;
Iwpl = (HMT~=0); Inwpl = (HMT==0);  % Two variables that define when workplace employment is positive or not positive (zero)
Irsd = (HRT~=0); Inrsd = (HRT==0);  % Two variables that define when residence employment is positive or not positive (zero)
nto = sum(Iwpl);                    % The number of observations with positive workplace employment
nfrom = sum(Irsd);                  % The number of observations with positive residence employment 
EHMT = HMT(Iwpl);                   % Generates a vector of workplace employment only containing observations with positive workplace employment (it has nto observations)
EHRT = HRT(Irsd);                   % Generates a vector of residence employment only containing observations with positive residence employment (it has nfrom observations)
Eomega = omega(Iwpl);               % Generates a vector of guesses of transformed wages for observations with positive workplace employment (it has nto observations)
Edistvar=distvar(Iwpl,Irsd);        % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom observations)
                                    % Notice that in this matrix, workplaces are rows and residences are columns 
clear distvar                       % clear variable to save memory

% Commuting costs;
Ecc=exp(kappaeps.*Edistvar);        % Generates an iceberg cost transport cost matrix based on the travel time matrix
clear Edistvar                      % clear variable to save memory
  
% Initializations;
gaptol=0;                           % This is the max difference between predicted and observe workplace employment that will be tolerated

% *************************;
% **** Calibrate Wages ****;
% *************************;

% xtic = tic();
display('>>>> Calibrating omegas <<<<');

% Staring OUTER loop; x is a counter; stopping rule is counter reaches 500
x=1;
while x<500; 

% Commuting probabilities;
    % Rows are workplace locations; we have nto rows
    % Columns are residence locations; we have nfrom columns
Eomegamat=repmat(Eomega,1,nfrom);   % we create a nto x nfrom matrix; each workplace employment (row, commuting destination) is assigned to all residence locations (commuting origins)
Ecnummat=Eomegamat./Ecc;            % We normalize destination transformed wages (at j) in commuting cost from each commuting origin (at i). Corresponds to computing the numerator in Eq. S.44 in the ARSW supplement
clear Eomegamat                     % clear variable to save memory
Ecdenom=sum(Ecnummat);              % Column-wise summation of commuting-cost deflated transformed wages. We sum over all workplaces for each residence. Corresponds to the the denominator in Eq. S.44 in the ARSW supplement. It has 1 x from elements, one for each residence
Ecdenommat=repmat(Ecdenom,nto,1);   % We repeat the Ecdenommat nto times, i.e. we create a matrix where each of the nto rows has the same elements in all nfrom colums
Ecprob=Ecnummat./Ecdenommat;        % We create the conditional commuting probabilities. Each element is the ratio of the commuting-cost deflated transformed wages over---the numerator in Eq. S.44---over the multilateral resitence term---the denominator in Eq. S.44. 
clear  Ecnummat Ecdenommat          % clear variable to save memory
% Check commuting probabilities sum to one within each coloumn;
test=sum(Ecprob);                   % Column-wise summation, returns sum of all commputing probability from one residence to all workplaces
mntest=mean(test);                  % The mean across all columns / residences => must be 1

% Employment;
EHMC=Ecprob*EHRT;                   % We perform a matrix multiplication of the nto x nfrom dimension commuting probability matrix and the nfrom x 1 vector of residence employment. The result is a nto x 1 vector of predicted workplace employment. This is the summation over i in eq. S.44.

% Define tolerance;
EHMT_r=round(EHMT.*100000);         % We scale up actual and predicted workplace employment by 1,000,000 before rounding
EHMC_r=round(EHMC.*100000);         % When we compute errors (gaps), our test will be more demanding (else it would be easier to have zero gaps)

% Check if have found the equilibrium omega vector;
gap=EHMC_r-EHMT_r;                  % Compute scaled gaps                 
gap=abs(gap);                       % Convert them into absolute gaps
gap=max(gap);                       % Keep the largest gap as one scalar
%[x mean(Eomega) gap]

% Use an if condition that evaluates the gap The stopping rule is that
% gap tolarance is reached (Set to zero in line 68)

% If gap tolerance is reached
    if gap<=gaptol;                    
        x=10000;                        % Set large x which will end the outer loop given the stopping rule in line 79
        display('>>>> Wage System Converged <<<<');
        wconverge=1;                    % Set convergence dummy to 1 to indicate convergence
    else 
% If fap tolerance is exceeded we need to update our guesses   
    Eomega_e=(EHMT./EHMC).*Eomega;      % Create new guesses of transformed wages. We inflate current guess by the ratio of observed workplace employment over predicted workplace employment
                                        % Intuitively, we will need greater transformed wages to rationalise data if we observe more workplace employment than we predict
    
    % Use an if condition to check if we obtain 'Not a Number' (NaN) values 
    if isnan(gap)==1;                   % If we get NaN, we need to fall back on random guesses to continue
    Eomega_e=0.95+(1.05-0.95).*rand(nto,1); % Generates a nto x 1 array of random updated guesses of transformed wages. The random values we raw range from 0 to one and have a mean of 0.5. Multplied by 0.1 we get a mean of 0.05. This is added to 0.95 to get the mean of one.
    Eomega=0.95+(1.05-0.95).*rand(nto,1);   % Generates a nto x 1 array of random old guesses of transformed wages using the same procedure.
    end;                                % This if-condition is a safetey net in case something goes wrong

    % If we do not have an NaN problem, we can proceed to update our
    % guesses to a weighted combination of old and new guesses
    Eomega=(0.5.*Eomega_e)+(0.5.*Eomega);   % Here we use a weight of 0.5

    % Choose units for omegas such that mean omega=1; Omega is identified
    % up to a constant since it is the numerator and the denominator.
    % Therefore any level we optain with the solver is arbitrary.
    mn_Eomega=mean(Eomega);             % We compute the mean of omega
    Eomega=Eomega./geomean(Eomega);     % We normalize omega by its geometric mean (the geometric mean is scale invariant)
    
end;    % End of if-else condition evaluating if gap tolerance is reached
x=x+1;  % Set counter to next value
gap
end;
%  End of OUTER loop. Outer loop ends after 500 iterations (see line 74 even if gap tolerance is violated)  


% Fill in non-zero entries; We have solved transformed wages for places
% with positive workplace employment only. But we have locations with zero
% workplace employment (if there is positive workplace or residence
% employment in any year)

omout=zeros(noj,1);                     % Create an N x 1 vector of zeros
omout(Iwpl)=Eomega;                     % Now we replace the values with the nto x 1 vector of transformed wages for locations where we have positive employment (Iwpl=1) 
                                        % Places with zero employment have a transformed wage of zero. This is theory-consistent since a zero-wage rationalizes zero employment
HMC=zeros(noj,1);                       % Create an N x 1 vector of zeros
HMC(Iwpl)=EHMC;                         % Now we replace the values with the nto x 1 vector of predicted workplace employment for locations where we have positive employment (Iwpl=1) 
                                          

% Fill in zeros in commuting probabilities;
cprob=zeros(noj,noj);                   % Create an N x N vector of zeros
cprob(Iwpl,Irsd)=Ecprob;                % Now replace the values with the nto x nfrom matrix of predicted workplace employment for the bilaterals with positive workplace and residence employment.

% Check we have an equilibrium;
%testcommute=HMT-(cprob*HRT);

    
    
    

