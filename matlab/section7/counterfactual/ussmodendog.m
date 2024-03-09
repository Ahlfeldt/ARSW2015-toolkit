%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB program file for the toolkit for             %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A verion of this program is part of the replication directory       %%%
%%% This version has been commented                                     %%%
%%% It has also been modified to be more memory efficient               %%%
%%% It has also been extended to record convergence paths               %%%
%%% The updating rule has been slightly altered to be more transparent  %%%
%%% and consistent with the structure of the model                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This function solves for endogenous model outcomes for given
%%% primitives, accounting for endogenous agglomeration effects
%%% It holds U?bar constant and corresponds to the openßcitz case
function [endog,ucprob,HH,Ubar,converge,convergencepath] = ussmodendog(param,fund,distvar,noj,UU)
        % This program uses the following inputs
        % param is a 1 x k vector of parameter values
        % fund is a 1 x k vector of N x 1 vectors of exogenous fundamentals
            % and starting values of endogenous variables
                % Fundamenal productivity, fundamental amenity, density of development, area,
                % floor space prices, workplace employment, residence
                % employment, wage, total income, commercial floor space share 
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n    
        % UU is the target reservation utility. Usually the one recovered
        % using ubar.m
    % This program produces the following outputs
        % endog is a vector of endogenous variables
            % Wage, total income, commercial floor space share, output,
            % residential rent, commercial rent, workplace employment,
            % residence employment
        % ucprob 
        % HH is total employment in the counterfactual
        % Ubar is expected utility in the counterfactual (held constant)
        % converge indicator that solver has converged
        % convergencepath Vector of gap measures by iteration
    % The names of the inputs need to correspond to objects that exist in
        % the workspace. The names of the outputs can be freely chosen
    % Basically, the function uses an iterative procedure starting from
    % guessed values of the following target variables: wages. floor space
    % prices and commercial floor space shares

% Clear large objects from previous runs to free up space
clear Cucprob CucprobCft ucprob06;                                      
    
% Declaring scalars as globals so they can be accessed and modified by any function    
global alpha beta kappa epsilon lambda delta eta rho;

% Start clock
xtic = tic();

% Extract paramer values from input
kappa=param(1);
epsilon=param(2);
lambda=param(3);
delta=param(4);
eta=param(5);
rho=param(6);

% Convergence indicator
converge=0;                                                                 % Set indicator to 0, will be set to 1 if convergence successful

% Gamma from Eq. (9)
gammaf=gamma((epsilon-1)./epsilon);

% Fundamentals;
% Column 1 : aBYw;
% Column 2 : bBYw; 
% Column 3 : VBYw;
% Column 4 : areaSYw;
% Column 5 : rentBYw; 
% Column 6 : empwplBYw; 
% Column 7 : emprsdBYw;
% Column 8 : wageBYw;
% Column 9 : vvBYw;
% Column 10 : thetaBYw;

% Read fundamentals and starting values that we feed into programme.
% Fundamentals have been inverted previously. We can feed observed values
% of endogenous variables for faster convergence, but the procedure also
% works with uniform guesses.
a=fund(:,1); b=fund(:,2); V=fund(:,3); K=fund(:,4);
QT=fund(:,5); HMT=fund(:,6); HRT=fund(:,7);
wT=fund(:,8); vvT=fund(:,9);
thetaT=fund(:,10);

% Use density of development and land area to compute total floor space
L=V.*(K.^0.75);                                                             % This uses Eq. (15)

% Only take dimensions with positives (positive workplace employment
% implies positive productivities. Same for residence employment and
% amenities
Ito = (a~=0); Into = (a==0);                                                % Two variables that define when productivity is positive or not positive (zero)
Ifrom = (b~=0); Infrom = (b==0);                                            % Two variables that define when amenity is positive or not positive (zero)
IcsA = ((a~=0)&(b==0));                                                     % Positive productivity and zero amenity; happens when blocks are fully specialized in production (positive workplace employment and zero residence employment)
IcsB = ((b~=0)&(a==0));                                                     % Positive amenity and zero productivity; happens when blocks are fully specialized in consumption (zero workplace employment and positive residence employment)
Iis = ((a~=0)&(b~=0));                                                      % Positive amenity and positive productivity, implying mixed land use
nto = sum(Ito);                                                             % Number of blocks with positive workplace employment
nfrom = sum(Ifrom);                                                         % Number of blocks with positive residence employment
nis = sum(Iis);                                                             % Number of blocks with positive residence and workplace employment
ea = a(Ito);                                                                % Vector of productivities, containing only blocks with positive workplace employment
eb = b(Ifrom);                                                              % Vector of amenities, containing only blocks with positive residence employment

% Initializations; creating guesses of endogenous objects marked by
% subscripts i. Predicted values will be marked by subscripts e.
Q_i=zeros(noj,1);                                                           % Create a J x 1 placeholder vector for guesses of residential floor space prices
Q_i(Ifrom)=QT(Ifrom);                                                       % Replace observations with guesses of residential floor space prices for locations with positive residence employment 
q_i=zeros(noj,1);                                                           % Create a J x 1 placeholder vector for guesses of commercial floor space prices
q_i(Ito)=QT(Ito);                                                           % Replace observations with guesses of commercial floor space prices for locations with positive workplace employment 
Q_e=Q_i;                                                                    % Set J x 1 value of predicted residential floor space prices for guesses. Will be updated later.
q_e=q_i;                                                                    % Set J x 1 value of predicted commercial floor space prices for guesses. Will be updated later.                                             
HM=HMT;                                                                     % New object containing observed workplace employment, will be updated while solving for counterfactual equilibrium
HR=HRT;                                                                     % New object containing observed residential employment, will be updated while solving for counterfactual equilibrium
HH=sum(HM);                                                                 % Total employment
wage_i=wT;                                                                  % Placeholder for guessed wages, contains observed wages to begin with 
wage_e=wage_i;                                                              % Placeholder for predicted wages
EHM=HM(Ito);                                                                % Workplace employment for locations with positive workplace emplyoment
EHR=HR(Ifrom);                                                              % Residence employment for locations with positive residence emplyoment 
EKM=K(Ito);                                                                 % Land area fro locations with positive workplace employment
EKR=K(Ifrom);                                                               % Land area fro locations with positive residence employment
Ewage_i=wage_i(Ito);                                                        % Vector of wage guesses for locations with postiive workplace employment
vv=vvT;                                                                     % New object containing total worker income, will be updated while solving for counterfactual equilibrium
theta_i=thetaT;                                                             % Set J x 1 value for guesses of floor space shares to observed shares. Will be updated later

% Commuting costs;
d_ij_eps=exp(-epsilon.*kappa.*distvar(Ifrom,Ito));                          % Component of route choice probabilities in Eq. (4)
d_ji_eps=d_ij_eps';                                                         % Transposed version of commuting cost matrix connecting residences to workplace (required for some calculations since sets of workplace and residence locations are not equal.

% Knowledge decay;
dd_ij=exp(-delta.*distvar(Ito,Ito));                                        % Spatial weights in productivity spillovers in Eq. (20)

% Residential decay;
cc_ij=exp(-rho.*distvar(Ifrom,Ifrom));                                      % Spatial weights in amenity spillovers in Eq. (21)

% Drop distvar to save space
clear distvar

% Endogenous amenities and productivities
EUps=dd_ij*(EHM./EKM);                                                      % Compute endogenos productivity Upsilon in Eq. (20) (we use matrix multiplication to execute the summation)
EA=ea.*(EUps.^lambda);                                                      % Compute total productivity using Eq. (20)
EOme=cc_ij*(EHR./EKR);                                                      % Compute endogenous amenity Omega in Eq. (21) (we use matrix multiplication to execute the summation)
EB=eb.*(EOme.^eta);                                                         % Compute total amenity using Eq. (21)

% CCompute output
Y=zeros(noj,1);                                                             % Placeholder for output, wil be updated later
Y(Ito)=EA.*(EHM.^alpha).*((theta_i(Ito).*L(Ito)).^(1-alpha));               % Compute actual output for locations with positive workplce employment using Eq. (10). Other blocks receive theory consistent zeros.    

% We generate a bunch placeholders for convergence measures fir 1000 iterations 
% that we can use to illustrate the covergence path
% Added by GA for didactic purposes
maxLDwagevec = zeros(1000,1);                                               % Max. log difference between guessed and predicted wages
maxLDqvec = zeros(1000,1);                                                  % Same for residential floor space prices
maxLDQvec = zeros(1000,1);                                                  % Same for commercial floor space prices                                                  
maxLDthetavec = zeros(1000,1);                                              % Same for commercial floorspace share
maxLDubarvec = zeros(1000,1);                                               % Same for U_bar
maxLDvec = zeros(1000,1);                                                   % Max. log diff across all locations and all outcomes
MSLEvec = zeros(1000,1);                                                    % Mean square log error acriss all locations and outcomes
Xvec = zeros(1000,1);                                                       % Identifier fot ith iteration

% Start the iterative procedure to solve for the equilibrium
% Stopping rule is when counter reaches 1000 (without convergence)
% Once guessess correspond to predicted values we set counter to 1000 (convergene)
% We need to find wages, floor space prices, commercial floor space shares
% and exogenous reservation utility
display('>>>> Solving the model <<<<');
x=1;
while x<1000;                                                               % Define stopping rule for outer loop; if counter reaches 1000. This can happen after 1000 iterations or because the if condition later on sets the counter to 1000 when convergence is reached

% We compute residence and workplace employment via commuting probabilities using guesses of wages and floor space prices and fundamental amenity;
% In the below matrices, we have living in i (rows) and working in j (cols);
% Notice that _i subscripts may indicate guesses (not residence locations, e.g. for wages)
EQQ=(Q_i(Ifrom).^(-(1-beta).*epsilon))*ones(1,nto);                         % Generate the first component of bilateral commuting probabilities in Eq. 4: floor space prices to the power of (1-beta) x -epsilon. Notice that the matrix multiplication by the 1 x nto vector of one enrues that we obtain a nfrom x nto matrix.
EBB=(EB*ones(1,nto)).^epsilon;                                              % Generates the next component of bilateral commuting probabilities: B to the power of epsilon. Again, the matrix multiplication by the 1 x nto vector of ones ensures that we obtain a nfrom x nto matrix. 
EWW=((Ewage_i*ones(1,nfrom))').^epsilon;                                    % Generates the next component of bilateral commuting probabilities, wage to the power of epsilon. Here, we start from a nto x 1 vector of wages. We use the matrix multiplication by the 1 x nfrom vector of ones to obtain a nto nfrom matrix, which we then transpose to have a nfrom x nto matrix.
Ephi_ij=d_ij_eps.*EBB.*EQQ.*EWW;                                            % Now we use the above components to compute the numerator of bilateral commuting probabilities in Eq. 4, phi_ij
Ephi_ji=Ephi_ij';                                                           % We transpose phi_ij to obtain a matrix that has workplaces in rows and residences in columns 
Ephi=sum(sum(Ephi_ij));                                                     % We sum over all bilaterals to obtain the denominator in Eq. 4. sum(Ephi_ij) retuns column sums, so we sum over the column sum to get the sum over all bilaterals
Ephi_i=sum(Ephi_ij'); Ephi_i=Ephi_i';                                       % We compute the numerators in Eq. 5: We transpose the phi_ij matrix so that residences are in columns. The row sum then returns sum(phi_ij) by i as a 1 x nrsd row vector. We transpose this vector to obtain a nrsd x 1 column vector.
Ephi_j=sum(Ephi_ij); Ephi_j=Ephi_j';                                        % We compute the numerators in Eq. 5: Summing over the columns of the nrsd x nwpl phi_ij matrix, returns sum(phi_ij) by j. We transpose the 1 x nwpl row vector to obtain a nwpl x 1 column vector.
Epp_ij=Ephi_ij./Ephi;                                                       % We divide the numerator by the denominator to obtain the commuting probabilities
Epp_ji=Ephi_ji./Ephi;                                                       % Transposed commuting probabilities with workplaces in rows and residences in columns    
Epp_i=Ephi_i./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a residence
Epp_j=Ephi_j./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a workplace
EHR=Epp_i.*HH;                                                              % Compute residence employment using predicted commuting probabilities for given guesses of wages and floor space prices
EHM=Epp_j.*HH;                                                              % Compute workplace employment predicted commuting probabilities for given guesses of wages and floor space prices
% Conditional commuting probabilities
Epp_iji=Ephi_ij./((Ephi_i*ones(1,nto)));                                    % To compute the conditional commuting probability, we normalize bilateral commuting probabilitiy by the residential choice probability. The multiplication by the 1 x nwpl row vector of ones serves to replicate the residential choice probabilities nwpl times so that we obtain an nrsd x nwpl matrix that can be used to normalize the bilateral commuting probability matrix
Epp_ijj=Ephi_ji./((Ephi_j*ones(1,nfrom)));                                  % We use the transposed matrix to compute the probability of living in i conditional on working in j
                                                                            % This matrix has commuting destinations in rows and commuting origins in columns, it is an nwpl x nres matrix
% Compute predicted reservation utility implied by Phi
Ubar=gammaf.*(Ephi.^(1./epsilon));                                          % Using Eq. (9)

% Knowledge spillovers;
EUps=dd_ij*(EHM./EKM);                                                      % Update Upsilon from Eq. (20) using current workpace employment
EA=ea.*(EUps.^lambda);                                                      % Update total productivity from Eq. (20) using current workpace employment

% Residential spillovers;
EOme=cc_ij*(EHR./EKR);                                                      % Update Omega from Eq. (21) using current workpace employment
EB=eb.*(EOme.^eta);                                                         % Update total amenity from Eq. (21) using current workpace employment

%%% Now we use guesses of commerical floor space shares and predicted
%%% employment (based on guesses of wages and floor space) to predict
%%% output
Y(Ito)=EA.*(EHM.^alpha).*((theta_i(Ito).*L(Ito)).^(1-alpha));               % Compute predicted output for given guesses of floor space prices and wages for locations with positive workplce employment using Eq. (10). Other blocks receive theory consistent zeros.

% Total worker income;
Evv=Epp_ijj'*(Ewage_i.*EHM);                                                % We compute expected working income according to Eq. S.20 using guesses of wages
vv(Ifrom)=Evv;                                                              % We write the total worker income for blocks with positive residence employment into the n x1 vector, those with resident employment obtain a theory-consistent zero total worker income

%test=(Epp_ij.*((Ewage_i*ones(1,nfrom))')).*HH;
%test=sum(test'); test=test';
%test=test-Evv;

%%% Now we use predicted output and workplace employment to predict wages %
Ewage_e=(alpha.*Y(Ito))./EHM;                                               % We use the first order condition for profit maximization with respect to labour for the production function in Eq. (10). This way, we obtain predicted wages, using predicted output (based on guessed theta and predicted emplyoment based on guessed wages) and workplace employment based on guesses of wages and floor space prices 

%%% Now we predict floor space prices using guesses of theta and predicted
%%% output that depends on guesses of theta and predicted wages that
%%% depends on guesses of wages and floor space prices

% Completely specialized commercial land prices;
q_e(IcsA)=((1-alpha).*Y(IcsA))./(theta_i(IcsA).*L(IcsA));                   % We use the first order-condition of profit-maximization, take the derivative of the profit function (Eq. 10 - minus costs for inputs) and maximize with respect to floor space prices q to get q as a function of output Y and commercial floor space (theta x LD)         

% Completely specialized residential land prices;
Q_e(IcsB)=((1-beta).*vv(IcsB))./((1-theta_i(IcsB)).*L(IcsB));               % We use the Marschallian demand function derived from Eq. (1). Expenditure on housing Q*LR = Q*(1-theta)*LD must equal the expenditure share on housing (1-beta) times total income (vv).

% Incompletely specialized land prices;
q_e(Iis)=(((1-alpha).*Y(Iis))+((1-beta).*vv(Iis)))./L(Iis);                 % For imperfectly specialized blocks, the rent is the weighted average, i.e. theta*q +(1-theta)*Q
Q_e(Iis)=(((1-alpha).*Y(Iis))+((1-beta).*vv(Iis)))./L(Iis);                 % For imperfectly specialized blocks, the rent is the weighted average, i.e. theta*q +(1-theta)*Q

%%% Now we predict commercial floor space shares using predicted values of
%%% rents and output

% Solve for theta for which the land market clears;
theta_e=theta_i;                                                            % Completely specialized blocks never change since amenity or productivity are zero, so predicted theta correponds guesses (initial values). We only update mixed use shares in the next line
theta_e(Iis)=((1-alpha).*Y(Iis))./(Q_e(Iis).*L(Iis));                       % We use the factor input demand function derived from Eq. (10). Commercial floor space input LM is (1-alpha)*Y / q. Theta is LM/LD. Combining the two we can predict the commerical floor space share for imperfectly specialized blocks using predicted outputs and rents

%%% Now we have closed to loop from guesses to predicted values of the 
%%% target variables wage, floor space price, floor space share and reservation utility (here impose a fixed level);
%%% We will update our guesses using predicted values until they converge.
%%% To this end, we need to define convergence criteria
wage_i(Ito)=Ewage_i;                                                        % Update N x 1 vector of guessed wages for locations with positive workplace employment (the others have zero anyways)
wage_e(Ito)=Ewage_e;                                                        % Update N x 1 vector of predicted wages for locations with positive workplace employment (the others have zero anyways)                                                       
wage_i_r=round(wage_i.*100);                                                % Round guessed wages to second digit - implicitly defines convergence precision as the stopping rule will be when guessed and predicted wages equate
wage_e_r=round(wage_e.*100);                                                % Round predicted wages to second digit - implicitly defines convergence precision as the stopping rule will be when guessed and predicted wages equate                                             
q_i_r=round(q_i.*100);                                                      % Similarly round floor space rents to second digit in this and the next three lines
q_e_r=round(q_e.*100);
Q_i_r=round(Q_i.*100);
Q_e_r=round(Q_e.*100);
theta_i_r=round(theta_i.*100);                                              % Similarly round commercial floor space share to second digit in this and the next line
theta_e_r=round(theta_e.*100);
UU_e=ones(noj,1).*Ubar;                                                     % Create a N x 1 vector of predicted reservation utility (will be useful in updating total employment below)
UU_i=ones(noj,1).*UU;                                                       % Create a N x 1 vector of imposed (as an input) reservation utility (will be useful in updating total employment below)
UU_e_r=round(UU_e.*100);                                                    % Round round predicted and imposed reservation utility to second digit in this and the next line
UU_i_r=round(UU_i.*100);

% Vector of differences between predicted and guessed values to illustrate
% convergence in real time; we will use a different vector
% [x; Ubar./UU; mean(abs(wage_i_r-wage_e_r)); mean(abs(q_i_r-q_e_r)); mean(abs(Q_i_r-Q_e_r)); mean(abs(theta_i_r-theta_e_r))]

%%% Now we evaluate if convergence is reached. 
%%% The stopping rule is that (rounded) guesses of target variables must
%%% correspond to rounded predicted values.
if (wage_e_r==wage_i_r) & (q_e_r==q_i_r) & (Q_e_r==Q_i_r) & (theta_e_r==theta_i_r) & (UU_e_r==UU_i_r);
    display('>>>> Convergence Achieved <<<<');
    display('Elapsed minutes');    
    xtic=toc(xtic);
    xtic./60
    x=10000;                                                                % If condition is (convergence) is fulfilled, we set the counter x to a large value. This executes the stopping rule in the outer loop
    converge=1;                                                             % We set the convergence indicator to one so that we know that we have converged
else;

%%% If we have not converged, we need to update our guesses %%%
%%% The remaining part of the code within the loop has been amended by GA
%%% to achieve faster convergence and save the convergence path

% We compute the mean log squared error to track the overall convergence of
% our target variables
MLSE=1000000./5.*(mean(abs(log(wage_i_r+1)-log(wage_e_r+1)).^2) + mean(abs(log(q_i_r+1)-log(q_e_r+1)).^2) + mean(abs(log(Q_i_r+1)-log(Q_e_r+1)).^2) + mean(abs(log(theta_i_r+1)-log(theta_e_r+1)).^2) + abs(log(Ubar)-log(UU)).^2 ) ;

% We compute the maximum log difference between the guessed and
% predicted value across observations in any of the targe variables
maxLDwage = max(abs(log(wage_i_r+1)-log(wage_e_r+1))) ;                     % wage
maxLDq = max(abs(log(q_i_r+1)-log(q_e_r+1)));                               % commercial floor space price 
maxLDQ = max(abs(log(Q_i_r+1)-log(Q_e_r+1)));                               % residential floor space price
maxLDtheta = max(abs(log(theta_i_r+1)-log(theta_e_r+1))) ;                  % commercial floor space share
maxLDUbar = abs(log(Ubar)-log(UU));                                         % reservation utility
maxLD = max([maxLDwage,maxLDq,maxLDQ,maxLDtheta,maxLDUbar]);                % The largest log difference across all target variables
[x;maxLDwage;maxLDq;maxLDQ;maxLDtheta;maxLDUbar;MLSE]                       % Output convergence measures to evaluate convergence in real time

% We read the convergence measures into vectors that we inspect later on
Xvec(x)=x;                                                                  % The iteration, the running variable
MSLEvec(x)=MLSE;                                                            % Mean log squared error
maxLDwagevec(x) = maxLDwage;                                                % Max log diff between guessed and predicted wages
maxLDqvec(x) = maxLDq;                                                      % Max log diff between guessed and predicted commercial floor space prices
maxLDQvec(x) = maxLDQ ;                                                     % Max log diff between guessed and predicted residential floor space prices
maxLDthetavec(x) = maxLDtheta;                                              % Max log diff between guessed and predicted commercial floor space shares
maxLDubarvec(x) = maxLDUbar;                                                % Max log diff between guessed and predicted reservation utility
maxLDvec(x) = maxLD;                                                        % Max log diff between guessed and predicted value across all target variables

%%% Now we update our guesses with subscript i using a weighted combination 
%%% of old guesses and predicted values (conditional on guesses) with subscript e

%ConvWeight_w = 0.25;
%ConvWeight_Q = 0.25;
%ConvWeight_theta = 0.25;
%ConvWeight_ubar = 0.25;

% We choose a weight depending on how far we are from convergence
if maxLDwage < 0.1                                                          % If we are far from convergence                                        
    ConvWeight_w = 0.25;                                                    % We give this weight to the predicted value
else ConvWeight_w = 0.5;                                                    % Else, we use this weight. It can be beneficial to choose a smaller weight if the the algorithm is bouncing too much close convergence
end
if maxLDQ < 0.1                                                             % If we are far from convergence                                        
    ConvWeight_Q = 0.25;                                                    % We give this weight to the predicted value
else ConvWeight_Q = 0.5;                                                    % Else, we use this weight. It can be beneficial to choose a smaller weight if the the algorithm is bouncing too much close convergence
end
if maxLDtheta < 0.1                                                         % If we are far from convergence                                        
    ConvWeight_theta = 0.25;                                                % We give this weight to the predicted value
else ConvWeight_theta = 0.5;                                                % Else, we use this weight. It can be beneficial to choose a smaller weight if the the algorithm is bouncing too much close convergence
end

% We update our guesses to weighted combinations of predicted values and
% old guesses
Ewage_i=(ConvWeight_w.*Ewage_e)+((1-ConvWeight_w).*Ewage_i);                % Wages
q_i=(ConvWeight_Q.*q_e)+((1-ConvWeight_Q).*q_i);                            % Floor commercial space prices
Q_i=(ConvWeight_Q.*Q_e)+((1-ConvWeight_Q).*Q_i);                            % Floor residential space prices
theta_i=(ConvWeight_theta.*theta_e)+((1-ConvWeight_theta).*theta_i);        % Commercial floor space share
%HH=(Ubar./UU).*HH;                                                         % This ad-hoc adjustment inspired by Eq. (9). We increase total employment if the utility in the city is higher than the user-specified target UU
HH_up=(Ubar./UU).^epsilon.*HH;                                              % In this variant of the updating rule we exploit the stucture of the model more explicitly. Linking Eq. (9) to the total employment equation defining H on p. 18 in the supplement, it is evident that utility in equilibrium must scale in total employment at an elasticity of (1/epsilon). By implication, employment must scale in total employment at an elasticity of epsilon which we exploit in this updating rule. 
HH=0.05*HH_up+0.95*HH                                                       % Since epsilon is relatively large, changes to predicted total employment can be large even if Ubar/UU is close to one. This can lead to bouncing and slow down convergence of the solver. Therefore, we assign a relatively small weights to the predicted to predicted total employment when updating total employment from iteration to iteration.
%HH=(0.25.*(Ubar./UU).*HH)+(0.75.*HH);

% Update loop;
% if (x==100) | (x==200) | (x==300) | (x==400) | (x==500) | (x==600) | (x==700) | (x==800) | (x==900) | (x==1000);
%     display('Iteration');
%     display(x);
%     display('Elapsed seconds');
%     display(toc(xtic));
% end;
x=x+1;
end;                                                                        % If-else condition for setting x to stopping rule ends
end;                                                                        % Outer loop ends

display('>>> Check Equilibrium Conditions <<<<');
test_j=EHM-(Epp_iji'*EHR);                                                  % Conditional choice probabilities correctly predict workplace employment
mean(test_j)
test_i=EHR-(Epp_ijj'*EHM);                                                  % Conditional choice probability correctly predict residence employment
mean(test_i)
test_i=EHR-(Epp_i.*HH);                                                     % Choice probabilities correctly predict residence employment
mean(test_i)
test_j=EHM-(Epp_j.*HH);                                                     % Choice probabilities correctly predict workplace employment
mean(test_j)

% Fill in non-zero entries;
HM=zeros(noj,1);                                                            % Placeholder for final workplace employment Nx1 vector
HM(Ito)=EHM;                                                                % Assign converged predicted values to locations with positive workplace employment, other remain theory-consistent zeros
HR=zeros(noj,1);                                                            % Placeholder for final residence employment Nx1 vector
HR(Ifrom)=EHR;                                                              % Assign converged predicted values to locations with positive residence employment, other remain theory-consistent zeros
A=zeros(noj,1);                                                             % Placeholder for total productivity Nx1 vector
A(Ito)=EA;                                                                  % Assign total productivity to locations with positive workplace employment, other remain theory-consistent zeros
B=zeros(noj,1);                                                             % Placeholder for total amenity Nx1 vector
B(Ifrom)=EB;                                                                % Assign total amenity to locations with positive residence employment, other remain theory-consistent zeros

% Fill in zeros commuting probabilities;
phi_ij=zeros(noj,noj);                                                      % N x N placeholder for unconditional commuting probabilities
phi_ij(Ifrom,Ito)=Ephi_ij;                                                  % Assign predicted commuting probabilities to routes with positive workplace and residence employment. Rest remains with theory-consistent zeros        
phi=sum(sum(phi_ij));                                                       % We sum over all bilaterals to obtain the denominator in Eq. 4. sum(Ephi_ij) retuns column sums, so we sum over the column sum to get the sum over all bilaterals
phi_i=sum(phi_ij'); phi_i=phi_i';                                           % We compute the numerators in Eq. 5: We transpose the phi_ij matrix so that residences are in columns. The row sum then returns sum(phi_ij) by i as a 1 x nrsd row vector. We transpose this vector to obtain a nrsd x 1 column vector.                                    
phi_j=sum(phi_ij); phi_j=phi_j';                                            % We compute the numerators in Eq. 5: Summing over the columns of the nrsd x nwpl phi_ij matrix, returns sum(phi_ij) by j. We transpose the 1 x nwpl row vector to obtain a nwpl x 1 column vector.
pp_i=phi_i./phi;                                                            % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a residence
pp_j=phi_j./phi;                                                            % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a workplace
phi_ji=phi_ij';                                                             % Transpose phi_ij matrix to connect workplaces to residences
pp_ij=phi_ij./phi;                                                          % Compute bilateral choice probabilities in Eq. (4)
pp_ji=phi_ji./phi;                                                          % Compute bilateral choice probabilities in Eq. (4)

temp=phi_i*ones(1,noj);                                                     % Cmompute conditional choice probabilities in Eq. (6)
pp_iji=phi_ij;
pp_iji(Ifrom,Ito)=pp_iji(Ifrom,Ito)./temp(Ifrom,Ito);
 
temp=phi_j*ones(1,noj);                                                     % Cmompute conditional choice probabilities conditional on workplace
pp_ijj=phi_ji;
pp_ijj(Ito,Ifrom)=pp_ijj(Ito,Ifrom)./temp(Ito,Ifrom);

display('>>> Check Equilibrium Conditions Again <<<<');
test_j=HM-(pp_iji'*HR);                                                     % Check that we correctly predict workplace employment form residence employment
mean(test_j)
test_i=HR-(pp_ijj'*HM);                                                     % Check that we correctly predict residence employment form worplace employment
mean(test_i)
test_i=HR-(pp_i.*HH);                                                       % Check that we correctly predict residence employment from residence choice choice probability
mean(test_i)
test_j=HM-(pp_j.*HH);                                                       % Check that we correctly predict workplace employment from workplace choice choice probability
mean(test_j)

% Unconditional commuting probabilities;
ucprob=zeros(noj,noj);                                                      % Placeholder for unconditional commuting probabilities
ucprob(Ifrom,Ito)=Epp_ij;                                                   % Write unconditional commuting probabilities for routes with positive commuting

% Crent;
Crent=zeros(noj,1);                                                         % Placeholder for rents
Crent(Ito)=q_i(Ito);                                                        % Write commercial rents for places with workplace employment
Crent(Ifrom)=Q_i(Ifrom);                                                    % Write residential rents for places with residence employment (all places have either workplace or residence employment)

% Endogenous variables and fundamentals output;
endog=[wage_i vv theta_i Y Q_i q_i HM HR Crent A B a b];

% Convergence path output
maxX = max(Xvec);
maxLDwagevec = maxLDwagevec(1:maxX);
maxLDqvec = maxLDqvec(1:maxX);
maxLDQvec = maxLDQvec(1:maxX); 
maxLDthetavec = maxLDthetavec(1:maxX);
maxLDvec = maxLDvec(1:maxX);
MSLEvec = MSLEvec(1:maxX);
Xvec = Xvec(1:maxX);
convergencepath = [maxLDwagevec maxLDqvec maxLDQvec maxLDthetavec maxLDvec MSLEvec Xvec];

display('<<<<<<<<<<< Equilibrium solved >>>>>>>>>>')



