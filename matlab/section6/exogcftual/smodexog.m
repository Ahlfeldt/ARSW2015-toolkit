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
%%% It has been modified to be more memory efficient                    %%%
%%% It has been extended to record convergence paths                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This function solves for endogenous model outcomes for given
%%% primitives. To this end, it uses several other functions. 
function [endog,ucprob,HH,converge,convergencepath,Ubar] = smodexog(param,fund,distvar,noj)
        % This program uses the following inputs
        % param is a 1 x k vector of parameters
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
        % endog is a vector of endogenous variables
            % ... 
        % ucprob 
        % converge indicator that solver has converged
        % convergencepath Vector of gap measures by iteration
    % The names of the inputs need to correspond to objects that exist in
        % the workspace. The names of the outputs can be freely chosen
    % Basically, the function uses an iterative procedure starting from
    % guessed values of the following target variables: wages. floor space
    % prices and commercial floor space shares


% Clear large objects from previous runs to free up space
clear Cucprob;                                  

% Declaring scalars as globals so they can be accessed and modified by any function
global alpha beta epsilon kappa kappaeps;

xtic = tic();                                                               % Start clock

% Extract paramer values from input
epsilon=param(1);
kappaeps=param(2);

% Gamma from Eq. 9
gammaf=gamma((epsilon-1)./epsilon);

% Indicator that suggests convergence. Switches to one if convergence
% successful
converge=0;

% Read vector of fundamentals and starting values of endogenous variables into individual vectors;
% Column 1 : ABYw;
% Column 2 : BBYw; 
% Column 3 : VBYw;
% Column 4 : areaSYw;
% Column 5 : rentBYw; 
% Column 6 : empwplBYw; 
% Column 7 : emprsdBYw;
% Column 8 : LM;
% Column 9 : LR;
% Column 10 : LD;
% Column 11 : wage;
% Column 12 : vv;

% Read fundamentals and starting values that we feed into programme.
% Fundamentals have been inverted previously. We can feed observed values
% of endogenous variables for faster convergence, but the procedure also
% works with uniform guesses.
A=fund(:,1); B=fund(:,2); V=fund(:,3); K=fund(:,4);
QT=fund(:,5); HMT=fund(:,6); HRT=fund(:,7);
LM=fund(:,8); LR=fund(:,9); LD=fund(:,10);
wage=fund(:,11); vv=fund(:,12);
HH=sum(HMT);

% Only take dimensions with positives (positive workplace employment
% implies positive productivities. Same for residence employment and
% amenities
Iwpl = (A~=0); Inwpl = (A==0);                                              % Two variables that define when productivity is positive or not positive (zero)
Irsd = (B~=0); Inrsd = (B==0);                                              % Two variables that define when amenity is positive or not positive (zero)
IcsA = ((A~=0)&(B==0));                                                     % Positive productivity and zero amenity; happens when blocks are fully specialized in production (positive workplace employment and zero residence employment)
IcsB = ((B~=0)&(A==0));                                                     % Positive amenity and zero productivity; happens when blocks are fully specialized in consumption (zero workplace employment and positive residence employment)
Iis = ((A~=0)&(B~=0));                                                      % Positive amenity and positive productivity, implying mixed land use
nwpl = sum(Iwpl);                                                           % Number of blocks with positive workplace employment
nrsd = sum(Irsd);                                                           % Number of blocks with positive residence employment    
nis = sum(Iis);                                                             % Number of blocks with positive residence and workplace employment
EA = A(Iwpl);                                                               % Vector of productivities, containing only blocks with positive workplace employment
EB = B(Irsd);                                                               % Vector of amenities, containing only blocks with positive residence employment

% Commuting costs;
Edistvar=distvar(Irsd,Iwpl);                                                % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom osbervations)
                                                                            % Notice that this time, we are having residence in rows and workplace in columsn, so we have an nfrom x nto matrix

d_ij_eps=exp(-kappaeps.*Edistvar);                                          % Compute exponential spatial weights declining in commuting decay
d_ji_eps=d_ij_eps';                                                         % Exponential spatial weights with workplace in rows and residence in columns (matters since original matrix has nfrom x nto dimensions)
Ecc=exp(-kappaeps.*distvar(Iwpl,Irsd));                                     % Exponential commuting weights from residence to workpalce
clear distvar                                                               % Drop variables 

% Initializations; creating guesses of endogenous objects marked by
% subscripts i. Predicted values are marked by subscripts e.
Q_i=zeros(noj,1);                                                           % Create a J x 1 placeholder vector for guesses of residential floor space prices
Q_i(Irsd)=QT(Irsd);                                                         % Replace observations with guesses of residential floor space prices for locations with positive residence employment 
q_i=zeros(noj,1);                                                           % Create a J x 1 placeholder vector for guesses of commercial floor space prices
q_i(Iwpl)=QT(Iwpl);                                                         % Replace observations with guesses of commercial floor space prices for locations with positive workplace employment 
Q_e=QT;                                                                     % Set J x 1 value of predicted residential floor space prices for guesses. Will be updated later.
q_e=QT;                                                                     % Set J x 1 value of predicted commercial floor space prices for guesses. Will be updated later. 
theta_i=LM./LD;                                                             % Use guesses of commercial and total floor space to set guesses for commercial floor space share
Ewage_i=wage(Iwpl);                                                         % Generate vector of guessed wages containing only locations with positive workplace employment
wage_i=zeros(noj,1);                                                        % Placeholder for guessed wages
wage_e=zeros(noj,1);                                                        % Placeholder for predicted wages
Y=zeros(noj,1);                                                             % Placeholder for output           
Y(Iwpl)=EA.*(HMT(Iwpl).^alpha).*(LM(Iwpl).^(1-alpha));                      % Compute actual output for locations with positive workplce employment using Eq. (10). Other blocks receive theory consistent zeros.    

display('>>>> Solving the model <<<<');
x=1;                                                                        % Set counter for iterations to 1
%MSLEvec=zeros(1000,1);
%vec=zeros(1000,1);
% We generate a bunch placeholders for convergence measures fir 1000 iterations 
% that we can use to illustrate the covergence path
% Added by GA for didactic purposes
maxLDwagevec = zeros(1000,1);                                               % Max. log difference between guessed and predicted wages
maxLDqvec = zeros(1000,1);                                                  % Same for residential floor space prices
maxLDQvec = zeros(1000,1);                                                  % Same for commercial floor space prices                                                  
maxLDthetavec = zeros(1000,1);                                              % Same for commercial floorspace share
maxLDvec = zeros(1000,1);                                                   % Max. log diff across all locations and all outcomes
MSLEvec = zeros(1000,1);                                                    % Mean square log error acriss all locations and outcomes
Xvec = zeros(1000,1);                                                       % Identifier fot ith iteration
% Back to original code

% Start the iterative procedure to solve for the equilibrium
% Stopping rule is when counter reaches 1000 (without convergence)
% Once guessess correspond to predicted values we set counter to 1000 (convergene)
% We need to find wages, floor space prices, and commercial floor space
% shares
while x<1000;                                                               % Define stopping rule for outer loop; if counter reaches 1000. This can happen after 1000 iterations or because the if condition later on sets the counter to 1000 when convergence is reached
    
    
% We compute residence and workplace employment via commuting probabilities using guesses of wages and floor space prices and fundamental amenity;
% In the below matrices, we have living in i (rows) and working in j (cols);
% Notice that _i subscripts may indicate guesses (not residence locations, e.g. for wages)
EQQ=repmat(QT(Irsd),1,nwpl);                                                % We assign guesses of residence floor space prices to all bilaterals. We take the vector of floor space prices, generate a new vector only for observations with positive residence employment, and replicates it nwpl times to generate a nrsd x nwpl matrix
EQQ=EQQ.^(-(1-beta).*epsilon);                                              % Generate the first component of bilateral commuting probabilities in Eq. 4: floor space prices to the power of (1-beta) x -epsilon 
EBB=repmat(EB,1,nwpl);                                                      % We assign residence amenities to all bilaterals. Takes the vector of amenities for observations with positive residence emplyoment and replicates it to generate an nrds x nwpl matrix
EBB=EBB.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities: B to the power of epsilon 
EWW=repmat(Ewage_i',nrsd,1);                                                % We assign guesses workplace wages to all bilaterals. We take the nwpl x 1 vector of adjusted wages, transpose it into an 1 x nwpl vector and replicate in nrsd rows to obtain an nrsd x nwpl matrix
EWW=EWW.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities, wage to the power of epsilon
Ephi_ij=d_ij_eps.*EBB.*EQQ.*EWW;                                            % Now we use the above components to compute the numerator of bilateral commuting probabilities in Eq. 4, phi_ij
Ephi_ji=Ephi_ij';                                                           % We transpose phi_ij to obtain a matrix that has workplaces in rows and residences in columns 
Ephi=sum(sum(Ephi_ij));                                                     % We sum over all bilaterals to obtain the denominator in Eq. 4. sum(Ephi_ij) retuns column sums, so we sum over the column sum to get the sum over all bilaterals
Ephi_i=sum(Ephi_ij'); Ephi_i=Ephi_i';                                       % We compute the numerators in Eq. 5: We transpose the phi_ij matrix so that residences are in columns. The row sum then returns sum(phi_ij) by i as a 1 x nrsd row vector. We transpose this vector to obtain a nrsd x 1 column vector.
Ephi_j=sum(Ephi_ij); Ephi_j=Ephi_j';                                        % We compute the numerators in Eq. 5: Summing over the columns of the nrsd x nwpl phi_ij matrix, returns sum(phi_ij) by j. We transpose the 1 x nwpl row vector to obtain a nwpl x 1 column vector.
Epp_ij=Ephi_ij./Ephi;                                                       % We divide the numerator by the denominator to obtain the commuting probabilities
Epp_ji=Ephi_ji./Ephi;                                                       % Transposed commuting probabilities with workplaces in rows and residences in columns           
Epp_i=Ephi_i./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a residence
Epp_j=Ephi_j./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a workplace
EHR=Epp_i.*HH;                                                              % Compute residence employment predicted commuting probabilities for given guesses of wages and floor space prices
EHM=Epp_j.*HH;                                                              % Compute workplace employment predicted commuting probabilities for given guesses of wages and floor space prices
Epp_iji=Ephi_ij./((Ephi_i*ones(1,nwpl)));                                   % To compute the conditional commuting probability, we normalize bilateral commuting probabilitiy by the residential choice probability. The multiplication by the 1 x nwpl row vector of ones serves to replicate the residential choice probabilities nwpl times so that we obtain an nrsd x nwpl matrix that can be used to normalize the bilateral commuting probability matrix
Epp_ijj=Ephi_ji./((Ephi_j*ones(1,nrsd)));                                   % We use the transposed matrix to compute the probability of living in i conditional on working in j
                                                                            % This matrix has commuting destinations in rows and commuting origins in columns, it is an nwpl x nres matrix
% Compute predicted reservation utility implied by Phi                                                                            
gammaf=gamma((epsilon-1)./epsilon);                                         % NEW LINE added for the teaching directory to return U_bar in the counterfactual
Ubar=gammaf.*(Ephi.^(1./epsilon));                                          % NEW LINE added for the teaching directory to return U_bar in the counterfactual
                                                                            
%%% Now we use guesses of commerical floor space shares and predicted
%%% employment (based on guesses of wages and floor space) to predict
%%% output
 Y(Iwpl)=EA.*(EHM.^alpha).*((theta_i(Iwpl).*LD(Iwpl)).^(1-alpha));         % Compute predicted output for given guesses of floor space prices and wages for locations with positive workplce employment using Eq. (10). Other blocks receive theory consistent zeros.

%%% Now we use predicted output and workplace employment to predict wages %
Ewage_e=(alpha.*Y(Iwpl))./EHM;                                              % We use the first order condition for profit maximization with respect to labour for the production function in Eq. (10). This way, we obtain predicted wages, using predicted output (based on guessed theta and predicted emplyoment based on guessed wages) and workplace employment based on guesses of wages and floor space prices 


% Total worker income;
% Evv=Epp_ijj'*(Ewage_i.*EHM);                                              % We transpose the matrix that has the probabilities of living in i conditional on working in j so that residences are in rows and workpalce are in columns. matrix multiplication of the nres x nwpl matrix by hte nres x 1 vector results in a nres x 1 vector
% GA alternative
Evv = Epp_iji * Ewage_i  ;                                                  % We compute expected working income according to Eq. S.20 using guesses of wages
Evv = Evv .* EHR ;                                                          % We multiply the expected worker income by the number of workers residing in that location (using guesses)
% GA alternative ends
vv(Irsd)=Evv;                                                               % We write the total worker income for blocks with positive residence employment into the n x1 vector, those with resident employment obtain a theory-consistent zero total worker income

% % Total worker income test;
% commute=Epp_ij.*HH;
% EWW=repmat(Ewage_i',nrsd,1);
% wcommute=commute.*EWW;
% test=zeros(noj,1);
% test(Irsd)=sum(wcommute,2);


%%% Now we predict floor space prices using guesses of theta and predicted
%%% output that depends on guesses of theta and predicted wages that
%%% depends on guesses of wages and floor space prices

% Completely specialized commercial block
q_e(IcsA)=((1-alpha).*Y(IcsA))./(theta_i(IcsA).*LD(IcsA));                  % We use the first order-condition of profit-maximization, take the derivative of the profit function (Eq. 10 - minus costs for inputs) and maximize with respect to floor space prices q to get q as a function of output Y and commercial floor space (theta x LD)         

% Completely specialized residential blocks;
Q_e(IcsB)=((1-beta).*vv(IcsB))./((1-theta_i(IcsB)).*LD(IcsB));              % We use the Marschallian demand function derived from Eq. (1). Expenditure on housing Q*LR = Q*(1-theta)*LD must equal the expenditure share on housing (1-beta) times total income (vv). 

% Incompletely specialized land prices;
q_e(Iis)=(((1-alpha).*Y(Iis))+((1-beta).*vv(Iis)))./LD(Iis);                % For imperfectly specialized blocks, the rent is the weighted average, i.e. theta*q +(1-theta)*Q
Q_e(Iis)=(((1-alpha).*Y(Iis))+((1-beta).*vv(Iis)))./LD(Iis);                % For imperfectly specialized blocks, the rent is the weighted average, i.e. theta*q +(1-theta)*Q

%%% Now we predict commercial floor space shares using predicted values of
%%% rents and output

% Solve for theta for which the land market clears;
theta_e=theta_i;                                                            % Completely specialized blocks never change since amenity or productivity are zero, so predicted theta correponds guesses (initial values). We only update mixed use shares in the next line
theta_e(Iis)=((1-alpha).*Y(Iis))./(q_e(Iis).*LD(Iis));                      % We use the facotr input demand function derived from Eq. (10). Commercial floor space input LM is (1-alpha)*Y / q. Theta is LM/LD. Combining the two we can predict the commerical floor space share for imperfectly specialized blocks using predicted outputs and rents

%%% Now we have closed to loop from guesses to predicted values of the 
%%% target variables wage, floor space price and floor space share;
%%% We will update our guesses using predicted values until they converge.
%%% To this end, we need to define convergence criteria
wage_i(Iwpl)=Ewage_i;                                                       % Update N x 1 vector of guessed wages for locations with positive workplace employment (the others have zero anyways)
wage_e(Iwpl)=Ewage_e;                                                       % Update N x 1 vector of predicted wages for locations with positive workplace employment (the others have zero anyways)                                                       
wage_i_r=round(wage_i.*100);                                                % Round guessed wages to second digit - implicitly defines convergence precision as the stopping rule will be when guessed and predicted wages equate
wage_e_r=round(wage_e.*100);                                                % Round predicted wages to second digit - implicitly defines convergence precision as the stopping rule will be when guessed and predicted wages equate                                             
q_i_r=round(q_i.*100);                                                      % Similarly round floor space rents to second digit in this and the next three lines
q_e_r=round(q_e.*100);
Q_i_r=round(Q_i.*100);
Q_e_r=round(Q_e.*100);
theta_i_r=round(theta_i.*100);                                              % Similarly round commercial floor space share to second digit in this and the next line
theta_e_r=round(theta_e.*100);

% Vector of differences between predicted and guessed values to illustrate
% convergence in real time; we will use a different vector
% [mean(abs(wage_i_r-wage_e_r)); mean(abs(q_i_r-q_e_r)); mean(abs(Q_i_r-Q_e_r)); mean(abs(theta_i_r-theta_e_r))]

%%% Now we evaluate if convergence is reached. 
%%% The stopping rule is that (rounded) guesses of target variables must
%%% correspond to rounded predicted values. 
if (wage_e_r==wage_i_r) & (q_e_r==q_i_r) & (Q_e_r==Q_i_r) & (theta_e_r==theta_i_r);
    display('>>>> Convergence Achieved <<<<');
    display('Elapsed seconds');    
    xtic=toc(xtic);
    xtic
    x=10000;                                                                % If condition is (convergence) is fulfilled, we set the counter x to a large value. This executes the stopping rule in the outer loop
    converge=1;                                                             % We set the convergence indicator to one so that we know that we have converged
else;
    
%%% If we have not converged, we need to update our guesses %%%
%%% The remaining part of the code within the loop has been amended by GA
%%% to achieve faster convergence and save the convergence path

x                                                                           % Output the counter to show the iteration of the loop
% We compute the mean log squared error to track the overall convergence of
% our target variables
MLSE=1000000./4.*(mean(abs(log(wage_i_r+1)-log(wage_e_r+1)).^2) + mean(abs(log(q_i_r+1)-log(q_e_r+1)).^2) + mean(abs(log(Q_i_r+1)-log(Q_e_r+1)).^2) + mean(abs(log(theta_i_r+1)-log(theta_e_r+1)).^2))

% We compute the maximum log difference between the guessed and
% predicted value across observations in any of the targe variables
maxLDwage = max(abs(log(wage_i_r+1)-log(wage_e_r+1))) ;                     % wage
maxLDq = max(abs(log(q_i_r+1)-log(q_e_r+1)));                               % commercial floor space price 
maxLDQ = max(abs(log(Q_i_r+1)-log(Q_e_r+1)));                               % residential floor space price
maxLDtheta = max(abs(log(theta_i_r+1)-log(theta_e_r+1))) ;                  % commercial floor space share
maxLD = max([maxLDwage,maxLDq,maxLDQ,maxLDtheta])                           % The largest log difference across all target variables

% We read the convergence measures into vectors that we inspect later on
Xvec(x)=x;                                                                  % The iteration, the running variable
MSLEvec(x)=MLSE;                                                            % Mean log squared error
maxLDwagevec(x) = maxLDwage;                                                % Max log diff between guessed and predicted wages
maxLDqvec(x) = maxLDq;                                                      % Max log diff between guessed and predicted commercial floor space prices
maxLDQvec(x) = maxLDQ ;                                                     % Max log diff between guessed and predicted residential floor space prices
maxLDthetavec(x) = maxLDtheta;                                              % Max log diff between guessed and predicted commercial floor space shares
maxLDvec(x) = maxLD;                                                        % Max log diff between guessed and predicted value across all target variables
    
%sum(abs(wage_i_r-wage_e_r).^2) + sum(abs(q_i_r-q_e_r).^2) + sum(abs(Q_i_r-Q_e_r).^2) + sum(abs(theta_i_r-theta_e_r).^2)

%%% Now we update our guesses with subscript i using a weighted combination 
%%% of old guesses and predicted values (conditional on guesses) with subscript e

% We choose a weight depending on how far we are from convergence
if maxLD < 0.5 %MLSE < 100                                                  % If we are far from convergence                                        
    ConvWeight = 0.5;                                                       % We give this weight to the predicted value
else ConvWeight = 0.5;                                                      % Else, we use this weight. It can be beneficial to choose a smaller weight if the the algorithm is bouncing too much close convergence
end

% We update our guesses
Ewage_i=(ConvWeight.*Ewage_e)+((1-ConvWeight).*Ewage_i);                    % Wages
q_i=(ConvWeight.*q_e)+((1-ConvWeight).*q_i);                                % Floor commercial space prices
Q_i=(ConvWeight.*Q_e)+((1-ConvWeight).*Q_i);                                % Floor residential space prices
theta_i=(ConvWeight.*theta_e)+((1-ConvWeight).*theta_i);                    % Commercial floor space share
%HH=Ephi;                                                                   % We want this to be closed-city model, so we do not update                     


x=x+1;
end;                                                                        % If-else condition for setting x to stopping rule ends
end;                                                                        % Outer loop ends

% display('>>> Check Equilibrium Conditions <<<<');
% test_j=EHM-(Epp_iji'*EHR);
% mean(test_j)
% test_i=EHR-(Epp_ijj'*EHM);
% mean(test_i)
% test_i=EHR-(Epp_i.*HH);
% mean(test_i)
% test_j=EHM-(Epp_j.*HH);
% mean(test_j)

% Fill in non-zero entries;
HM=zeros(noj,1);                                                            % Placeholder for final workplace employment Nx1 vector
HM(Iwpl)=EHM;                                                               % Assign converged predicted values to locations with positive workplace employment, other remain theory-consistent zeros
HR=zeros(noj,1);                                                            % Placeholder for final residence employment Nx1 vector
HR(Irsd)=EHR;                                                               % Assign converged predicted values to locations with positive residence employment, other remain theory-consistent zeros

% Fill in zeros commuting probabilities;
phi_ij=zeros(noj,noj);                                                      % N x N placeholder for unconditional commuting probabilities
phi_ij(Irsd,Iwpl)=Ephi_ij;                                                  % Assign predicted commuting probabilities to routes with positive workplace and residence employment. Rest remains with theory-consistent zeros        
phi=sum(sum(phi_ij));  
ucprob = phi_ij ./ phi;

% Crent;
Crent=zeros(noj,1);
Crent(Iwpl)=q_i(Iwpl);
Crent(Irsd)=Q_i(Irsd);

% Endogenous variables;
endog=[wage_i vv theta_i Y Q_i q_i HM HR Crent];

%%% GA code starts %%%
maxX = max(Xvec);
maxLDwagevec = maxLDwagevec(1:maxX);
maxLDqvec = maxLDqvec(1:maxX);
maxLDQvec = maxLDQvec(1:maxX); 
maxLDthetavec = maxLDthetavec(1:maxX);
maxLDvec = maxLDvec(1:maxX);
MSLEvec = MSLEvec(1:maxX);
Xvec = Xvec(1:maxX);
convergencepath = [maxLDwagevec maxLDqvec maxLDQvec maxLDthetavec maxLDvec MSLEvec Xvec];
%%% GA code ends%%%

