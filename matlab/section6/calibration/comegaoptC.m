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
%%% It has been modified to be more memory efficient                    %%%
%%% This version has been commented                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This function solves for transformed wages (omega) for given values of 
%%% workplace employment, residence employment, and bilateral commuting 
%%% costs
%%% Then, it uses solved transformed wages to recover adjusted productivity 

function [wgout,Aout,cprob,wconverge,HMC,gap] = comegaoptC(obsdata,distvar,noj,wgin);
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
        % wgin is vector of guesses of ADJUSTED wages (notice, comegaoptO
            % delivers TRANSFORMED wages
    % This program produces the following outputs
        % wgout is the solved equilibrium vector of ADJUSTED wages 
        % Aout is is a vector of inverted adjusted productivities
        % cprob is the matrix of condional commuting probabilities (rows
            % are workplaces, columns are residences)
        % wconverge is a dummy indicating that the solver has convered
        % HMC is predicted workplace employment (should the same as
            % observed workplace employment in case of perfect convergence)
        % gap is the largest difference between any pair of workplace
            % employment and predicted workplace employment
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen
        
    % Notice that unlike comegaoptO.m this function needs epsilon to be set
    % or estimated and saved as a global!
    
% Declaring scalars as globals so they can be accessed and modified by any function
global alpha beta epsilon kappa;

% Setting the convergence indicator to zero (if convergence is reached it
% will be replaced with 1
wconverge=0;

% Extracting variables from the four key input variables
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area

% Generate TRANSFORMED wages using guesses of ADJUSTED wages and epsilon
omega=wgin.^epsilon;

% Only take dimensions with positive employment;
Iwpl = (HMT~=0); Inwpl = (HMT==0);                                          % Two variables that define when workplace employment is positive or not positive (zero)
Irsd = (HRT~=0); Inrsd = (HRT==0);                                          % Two variables that define when residence employment is positive or not positive (zero)
nto = sum(Iwpl);                                                            % The number of observations with positive workplace employment
nfrom = sum(Irsd);                                                          % The number of observations with positive residence employment
EHMT = HMT(Iwpl);                                                           % Generates a vector of workplace employment only containing observations with positive workplace emplyoment (it has nto observations)
EHRT = HRT(Irsd);                                                           % Generates a vector of residence employment only containing observations with positive residence emplyoment (it has nfrom observations)
Eomega = omega(Iwpl);                                                       % Generates a vector of guesses of transformed wages for observations with positive workplace employment (it has nto observations)
Edistvar=distvar(Iwpl,Irsd);                                                % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom osbervations)
                                                                            % Notice that in this matrix, workplaces are rows and residences are columns 
clear distvar                                                               % clear variable to save memory

% Commuting costs;
Ecc=exp(epsilon.*kappa.*Edistvar);                                          % Generates an iceberg cost transport cost matrix based on the travel time matrix
clear Edistvar                                                              % clear variable to save memory

% Initializations;
gaptol=0;                                                                   % This is the max difference between predicted and observe workplace employment that will be tolerated

% *************************;
% **** Calibrate Wages ****;
% *************************;

% xtic = tic();
display('>>>> Calibrating omegas <<<<');

% Staring OUTER loop; x is a counter; stopping rule is counter reaches 500
x=1;
while x<500;

% The below lines have been rewritten to be more memory efficent. The 
% recovered wages and amenities are not exactly identical, therefore this
% alternative codes should only be used as a last resort
%%% New code by GA starts %%%
%{
clear Ecprob % we drop residuals from previous iterations in the loop to create space
% Assuming nto is the number of workplace locations and nfrom is the number of residence locations
% Assuming Eomega is a vector of length nto and Ecc is either a scalar or a vector of length nto

% Initialize Ecprob with zeros
Ecprob = zeros(nto, nfrom);

% Process each residence location
for col = 1:nfrom
    % Calculate commuting probabilities for each workplace location
    for row = 1:nto
        if isscalar(Ecc)
            Ecprob(row, col) = Eomega(row) / Ecc;  % If Ecc is a scalar
        else
            Ecprob(row, col) = Eomega(row) / Ecc(row);  % If Ecc is a vector
        end
    end
end

% Normalize the columns of Ecprob so that probabilities sum to one for each column
for col = 1:nfrom
    Ecprob(:, col) = Ecprob(:, col) / sum(Ecprob(:, col));
end

% Check if the columns of Ecprob sum to one
columnSums = sum(Ecprob);
% Now columnSums should be a row vector of ones if everything is correct
%%% New code by GA ends %%%
%}

% We use the original code  
% Commuting probabilities;
    % Rows are workplace locations; we have nto rows
    % Columns are residence locations; we have nfrom columns
Eomegamat=repmat(Eomega,1,nfrom);                                           % we create an nto x nfrom matrix; each workplace employment (row, commuting destination) is assigned to all residence locations (commuting origins)
Ecnummat=Eomegamat./Ecc;                                                    % We normalize destination transformed wages (at j) in commuting cost from each commuting origin (at i). Corresponds to computing the numerator in Eq. S.44 in the ARSW supplement
clear Eomegamat                                                              % clear variable to save memory
Ecdenom=sum(Ecnummat);                                                      % Column-wise summation of commuting-cost deflated transformed wages. We sum over all workplaces for each residence. Corresponds to the the denominator in Eq. S.44 in the ARSW supplement. It has 1 x from elements, one for each residence
Ecdenommat=repmat(Ecdenom,nto,1);                                           % We repeat the Ecdenommat nto times, i.e. we create a matrix where each of the nto rows has the same lements in all nfrom colums
Ecprob=Ecnummat./Ecdenommat;                                                % We create the conditional commuting probabilities. Each element is the ratio of the commuting-cost deflated transformed wages over---the numerator in Eq. S.44---over the multilateral resitence term---the denominator in Eq. S.44. 
clear  Ecnummat Ecdenommat                                                  % clear variable to save memory
    


% Check commuting probabilities sum to one in the columns of the;
% commuting matrix for each block of residence;
test=sum(Ecprob);                                                           % Column-wise summation, returns sum of all commputing probability from one residence to all workplaces
mntest=mean(test);                                                          % The mean across all columns / residences => must be 1
% Employment;
EHMC=Ecprob*EHRT;                                                           % We perform a matrix multiplication of the nto x nfrom dimension commuting probability matrix and the nfrom x 1 vector of residence employment. The result is a nto x 1 vector of predicted workplace employment. This is the summation over i in eq. S.44.
% Define tolerance;
EHMT_r=round(EHMT.*10000);                                                  % We scale up actual and predicted workplace employment by 1,000,000 before rounding
EHMC_r=round(EHMC.*10000);                                                  % When we compute errors (gaps), our test will be more demanding (else it would be easier to have zero gap)

% Check if have found the equilibrium omega vector;
gap=EHMC_r-EHMT_r;                                                          % Compute scaled gaps        
gap=abs(gap);                                                               % Convert them into absolute gaps
gap=max(gap)                                                                % Keep the largest gap as one scalar
%[x mean(Eomega) gap]

% Use an if condition that evaluates the gap. The stopping rule is that
% gap tolarance is reached (Set to zero in line 74)

% If gap tolerance is reached
if gap<=gaptol;
    x=10000;                                                                % Set large x which will end the outer loop given the stopping rule in line 85
    display('>>>> Wage System Converged <<<<');
    wconverge=1;                                                            % Set convergance dummy to 1 to indicate convergence
else 
 % If fap tolerance is exceeded we need to update our gesses      
    Eomega_e=(EHMT./EHMC).*Eomega;                                          % Create new gesses of transformed wages. We inflate current guess by the ratio of observed workplace emloyment over predicted workplace employment
                                                                            % Intuitively, we will need greater transformed wages to rationalize data if we observe more workplace employment than we predict
    
    % Use an if condition to check if we obtain 'Not a Number' (NaN) values 
    if isnan(gap)==1;                                                       % If we get NaN, we need to fall back on random guesses to continue
    Eomega_e=0.95+(1.05-0.95).*rand(nto,1);                                 % Generates a nto x 1 array of random updated guesses of transformed wages. The random values we raw range from 0 to one and have a mean of 0.5. Multplied by 0.1 we get a mean of 0.05. This is added to 0.95 to get the mean of one.                            
    Eomega=0.95+(1.05-0.95).*rand(nto,1);                                   % Generates a nto x 1 array of random old guesses of transformed wages using the same procedure.                               
    end;                                                                    % This if condition is a saftey net in case something goes wrong

    % If we of not have an NaN problem, we can proceed to update our
    % guesses to a weighted combination of old and new guesses
    Eomega=(0.5.*Eomega_e)+(0.5.*Eomega);

    % Choose units for omegas such that mean omega=1; Omega is identified
    % up to a constant since it the numerator and the denominator. Thefore
    % any level we uptain with the solver is arbitrary.
    Eomega=Eomega./geomean(Eomega);
    
end;                                                                        % End of if-else condition evaluating if gap tolerance is reached
x=x+1;                                                                      % Set counter to next value
end;
%  End of OUTER loop. Outer loop ends after 500 iterations (see line 74 even if gap tolerance is violated)  
    
% Fill in non-zero entries; We have solved transformed wages for places
% with positive workplace employment ony. But we have locations with zeor
% workplace employment (if there is positive workplace or residence
% employment in any year)

omega=zeros(noj,1);                                                         % Create an N x 1 vector of zeros
omega(Iwpl)=Eomega;                                                         % Now we replace the values with the nto x 1 vector of transformed wages for locations where we have positive employment (Iwpl=1) 
                                                                            % Places with zero employment have a tranformed wage of zero. This is theory-consistent since a zero-wage rationalizes zero employment
HMC=zeros(noj,1);                                                           % Create an N x 1 vector of zeros
HMC(Iwpl)=EHMC;                                                             % Now we replace the values with the nto x 1 vector of predicted workplace employment for locations where we have positive employment (Iwpl=1) 

% Fill in zeros commuting probabilities;
% The next two lines of original code will be replaced to be more memory
% efficient
%cprob=zeros(noj,noj);                                                       % Create an N x N vector of zeros
%cprob(Iwpl,Irsd)=Ecprob;                                                    % Now replace the values with the nto x nfrom matrix of predicted workplace employment for the bilaterals with positive workplace and residence employment.

%%% New code by GA start %%%
% Initialize cprob as a sparse matrix
cprob = sparse(noj, noj);

% Add values incrementally
for idx = 1:numel(Iwpl)
    cprob(Iwpl(idx), Irsd(idx)) = Ecprob(idx);
end
%%% New code by GA ends %%%

% Recover wages;
wgout=omega.^(1./epsilon);                                                  % Go back from TRANSFORMED wages to ADJUSTED wages for which are solving here
wgout(wgout>0)=wgout(wgout>0)./geomean(wgout(wgout>0));                     % Normalize positive wages by their geometric mean (zero wages represent theory-consistently zero productivity

% Recover productivity;
Aout=zeros(noj,1);                                                          % Generate 
Aout(Iwpl)=((QT(Iwpl)./(1-alpha)).^(1-alpha)).*((wgout(Iwpl)./alpha).^alpha); % We use Eq. S.48 (bottom one) solved for A_tilde to obtain this equation

% Check we have an equilibrium;
% testcommute=HMT-(cprob*HRT)

    
    
    

