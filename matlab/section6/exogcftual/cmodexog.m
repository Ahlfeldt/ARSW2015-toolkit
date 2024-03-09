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
%%% It has also been modified to be more memory efficient               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This function solves for transformed wages (omega) for given values of 
%%% workplace employment, residence employment, and bilateral commuting 
%%% costs
%%% Then, it uses solved transformed wages to recover adjusted productivity 

function [Aout,Bout,wout,ucprob,vv,HMC,HRC,CMA,Ephi,HH,ABconverge,mAgap,mBgap] = cmodexog(obsdata,distvar,noj,A,B)
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
        % A contains guesses of producitivities
        % B contains guesses of amenities
    % This program produces the following outputs
        % Aout is is a vector of inverted adjusted productivities
        % Bout is is a vector of inverted adjusted amenities
        % ucprob is the matrix of uncondional commuting probabilities 
            % residence in rows, workplace in columns
        % vv is total worker income
        % HMC is predicted workplace employment (should the same as
            % observed workplace employment in case of perfect convergence)
        % HRC is predicted residence employment (should the same as
            % observed residence employment in case of perfect convergence)
        % CMA is (residential) commuting market access
        % EPhi is the denominator in Eq. 4
        % HH is the sum of employment
        % ABconverge is an indicator that is one if there was convergence
        % mAgap and mBgap are convergence measures, should be zero when
            % there was convergence
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen
    
    % Intuitively, the programme starts from guesses of wages, amenities, and
    % productivities and updates them iteratively until an internally
    % consistent solutions is found. The program also adjusts amenities to
    % ensure that the population in the model matches the population in the
    % data
   
% Declaring scalars as globals so they can be accessed and modified by any function    
global alpha beta kappa epsilon;

% Setting the convergence indicator to zero (if convergence is reached it
% will be replaced with 1
ABconverge=0;

% Extracting variables from the four key input variables
    QT=obsdata(:,1);    % floor space prices
    HMT=obsdata(:,2);   % workplace employment
    HRT=obsdata(:,3);   % residence employment
    K=obsdata(:,4);     % geographic area

% gamma function;
gammaf=gamma((epsilon-1)./epsilon);                                         % We create the scalar gamma in the expected utiltiy Eq. (9)

% Only take dimensions with positive employment
Iwpl = (HMT~=0); Inwpl = (HMT==0);                                          % Two variables that define when workplace employment is positive or not positive (zero)
Irsd = (HRT~=0); Inrsd = (HRT==0);                                          % Two variables that define when residence employment is positive or not positive (zero)
nwpl = sum(Iwpl);                                                           % The number of observations with positive workplace employment
nrsd = sum(Irsd);                                                           % The number of observations with positive residence employment 
EHMT = HMT(Iwpl);                                                           % Generates a vector of workplace employment only containing observations with positive workplace emplyoment (it has nto observations)
EHRT = HRT(Irsd);                                                           % Generates a vector of residence employment only containing observations with positive residence emplyoment (it has nfrom observations)
wage = zeros(noj,1);                                                        % Placeholder for wage vector to be generated 
Ewage = wage(Iwpl);                                                         % Generates a vector of guesses of wages only containing observations with positive workplace employment
EA = A(Iwpl);                                                               % Generates a vector of adjusted productivities only containing obserations with positive workplace employment
EB = B(Irsd);                                                               % Generates a vector of adjusted amenities only containing obserations with positive residence employment
EKM = K(Iwpl);                                                              % Generates a vector of geographic area only only containing observations with positive workplace employment
EKR = K(Irsd);                                                              % Generates a vector of geographic area only only containing observations with positive residence employment    
HH=sum(EHMT);                                                               % Compute total workplace employment

% Commuting costs;
Edistvar=distvar(Irsd,Iwpl);                                                % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom osbervations)
clear distvar                                                               % Notice that this time, we are having residence in rows and workplace in columsn, so we have an nfrom x nto matrix
d_ij_eps=exp(-epsilon.*kappa.*Edistvar);                                    % Generates a spatial weight matrix with weights that exponentially falls in bilateral travel time between i and j at rate epsilon*kappa
clear Edistvar                                                              % clear variable to save memory 

% Initializations;
gaptol=0;                                                                    % This is the max difference between guessed and predicted values that will be tolerated

% *************************;
% **** Calibrate Wages ****;
% *************************;

display('>>>> Calibrating productivity and amenities <<<<');

% Staring OUTER loop; x is a counter; stopping rule is counter reaches 200
x=1;
while x<200;

%%% We begin by computing commuting probabilities (and hence workplace and
%%% residence employment) that is consistent with guesses of productivities
%%% and amenities and the resulting wages

% The following line has been brought up by GA so that wages consistent with
% current guesses of A are input into the commuting probabilities
Ewage=(((1-alpha)./QT(Iwpl)).^((1-alpha)./alpha)).*alpha.*(EA.^(1./alpha)); % We use Equation (12) which combines first-order condition and zero-profit conditions.
 
% The following lines represent new code by GA to compute the numerator of
% Eq. 4 more memory efficiently. It can be used in cases of severe memory
% shortage. However, it is significantly slower
%{
Ephi_ij = zeros(nrsd, nwpl);
% Compute Ephi_ij element-wise
EQT=QT(Irsd);
for i = 1:nrsd
    for j = 1:nwpl
        % Assuming QT(Irsd(i)), EB(i), Ewage(j), and d_ij_eps(i, j) are scalar values
        EQQ_element = EQT(i)^(-(1-beta) * epsilon);
        EBB_element = EB(i)^epsilon;
        EWW_element = Ewage(j)^epsilon;
        Ephi_ij(i, j) = d_ij_eps(i, j) * EQQ_element * EBB_element * EWW_element;  % Routes will be attractive if they offer high wages, high amenities and low rents, therefore attracting a larger share of commuters
    end
end
% New GA code ends
%}
% The below lines should be outcommented if the above lines are used
% Employments;
% Living in i (rows) working in j (cols);
EQQ=repmat(QT(Irsd),1,nwpl);                                                % We assign residence floor space prices to all bilaterals. We take the vector of floor space prices, generate a new vector only for observations with positive residence employment, and replicates it nwpl times to generate a nrsd x nwpl matrix
EQQ=EQQ.^(-(1-beta).*epsilon);                                              % Generate the first component of bilateral commuting probabilities in Eq. 4: floor space prices to the power of (1-beta) x -epsilon 
EBB=repmat(EB,1,nwpl);                                                      % We assign residence amenities to all bilaterals. Takes the vector of amenities for observations with positive residence emplyoment and replicates it to generate an nrds x nwpl matrix
EBB=EBB.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities: B to the power of epsilon 
EWW=repmat(Ewage',nrsd,1);                                                  % We assign workplace wages to all bilaterals. We take the nwpl x 1 vector of adjusted wages, transpose it into an 1 x nwpl vector and replicate in nrsd rows to obtain an nrsd x nwpl matrix
EWW=EWW.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities, wage to the power of epsilon
Ephi_ij=d_ij_eps.*EBB.*EQQ.*EWW;                                            % Now we use the above components to compute the numerator of bilateral commuting probabilities in Eq. 4, phi_ij
clear EQQ EWW EBB                                                           % Clear objects to save memory
% Outcomment until here

% Original code continues
Ephi=sum(sum(Ephi_ij));                                                     % We sum over all bilaterals to obtain the denominator in Eq. 4. sum(Ephi_ij) retuns column sums, so we sum over the column sum to get the sum over all bilaterals
Ephi_i=sum(Ephi_ij'); Ephi_i=Ephi_i';                                       % We compute the numerators in Eq. 5: We transpose the phi_ij matrix so that residences are in columns. The row sum then returns sum(phi_ij) by i as a 1 x nrsd row vector. We transpose this vector to obtain a nrsd x 1 column vector.
Ephi_j=sum(Ephi_ij); Ephi_j=Ephi_j';                                        % We compute the numerators in Eq. 5: Summing over the columns of the nrsd x nwpl phi_ij matrix, returns sum(phi_ij) by j. We transpose the 1 x nwpl row vector to obtain a nwpl x 1 column vector.
Epp_ij=Ephi_ij./Ephi;                                                       % We divide the numerator by the denominator to obtain the commuting probabilities                                                  
Epp_i=Ephi_i./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a residence
Epp_j=Ephi_j./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a workplace
EHRC=Epp_i.*HH;                                                             % Compute residence employment predicted commuting probabilities (should be equal to observed reidence employment)
EHMC=Epp_j.*HH;                                                             % Compute workplace employment predicted commuting probabilities (should be equal to observed reidence employment)


%%% Now define our convergence criterion. We want predicted workplace and
%%% residence employment to exactly match the values obeserved in data
Agap=EHMT-EHMC;                                                             % Difference between observed and predicted workplace employment
Bgap=EHRT-EHRC;                                                             % Difference between observed and predicted residence employment
Agap_r=round(Agap.*10000);                                                  % We round to the fifth digit, implicitly defines precision
Bgap_r=round(Bgap.*10000);                                                  % We round to the fifth digit, implicitly defines precision
mAgap=max(abs(Agap_r))                                                      % Compute absolute value of gap
mBgap=max(abs(Bgap_r))                                                      % Compute absolute value of gap
% [x mAgap mBgap];

%%% We stop the iterative procedure if our convergence criteria are met
if (mAgap<=gaptol) && (mBgap<=gaptol);
    display('>>>> Calibration Convergence Achieved <<<<');
    x=1000000;                                                              % Setting x to a large value executes the stopping rule when x is evaluated in the outer loop
    ABconverge=1;                                                           % Set the convergence indicator to 1 = success
else;

    %%% If there is no convergence, we need to update our guesses of amenities,
    %%% productivities, and, consequentially, wages. Intuitively. we increase
    %%% productivities where we underpredict workplace employment and amenities
    %%% where we underpredict residence employment
    EA_e=((EHMT./EHMC).^(1./epsilon)).*EA;                                      % Generate updated productivities: If ratio of actual over predicted workplace employment is > 1, we increase productivity, which will lead to larger predicted employment next round
    EB_e=((EHRT./EHRC).^(1./epsilon)).*EB;                                      % Generate updated amenities: If ratio of actual over predicted residence employment is > 1, we increase amenity, which will lead to larger residence employment next round

    % Use an if condition to check if we obtain 'Not a Number' (NaN) values 
    if (isnan(mBgap)==1);                                                       % If we get NaN, we need to fall back on random guesses to continue
    EA_e=0.95+(1.05-0.95).*rand(nwpl,1);                                        % Generates a nwpl x 1 array of random updated guesses of productivities. The random values we draw range from 0 to one and have a mean of 0.5. Multplied by 0.1 we get a mean of 0.05. This is added to 0.95 to get the mean of one.                              
    EB_e=0.95+(1.05-0.95).*rand(nrsd,1);                                        % Similarly generates a nrsd x 1 array of updated guesses of amenities
    Ewage=(((1-alpha)./QT(Iwpl)).^((1-alpha)./alpha)).*alpha.*(EA_e.^(1./alpha)); % Update wage to be consistent with the fall-back productivities
    EA=EA_e;
    EB=EB_e;
    end;

    %%% We update productivities and amenities to a weighted combination of
    %%% previous and updated guesses. If we fall back on the fall-back
    %%% fundamentals, the next two lines are inconsequential since EA =
    %%% EA_e and EB = AB_e
    EA=(0.5.*EA_e)+(0.5.*EA);                                                   % Update productivities
    EB=(0.5.*EB_e)+(0.5.*EB);                                                   % Update amenities

    % Choose units mean productivity =1; Productivities are identified
    % up to a constant since it the numerator and the denominator. Thefore
    % any level we uptain with the solver is arbitrary.
    EA=EA./geomean(EA);
    % Adjust level of amenities so that it is consistent with observed
    % population. We measure utility in units so that
    % (U_bar/gamma)^epsilon/H = 1. See supplement p. 17. This definition
    % ensures that Phi = H. To see this, compare H on supplement, p. 18 to
    % Phi, the denominator in Eq. (4). Thus, if the population in the data
    % (H) is greater than the population in the model (Phi), we increase
    % the amenity level to make the city more attractive and attract more
    % residents.
    EB=((HH./Ephi).^(1./epsilon)).*EB;                                      % If the ratio of the total population in the data over the population in the model is positive, we inflate amenities
x=x+1;
end;
end;

% Fill in non-zero entries;
Aout=zeros(noj,1);                                                          % Placeholder for N x 1 productivity vector
Aout(Iwpl)=EA;                                                              % Fill locations with positive workplace employment. Those without receive theory-consistent zeros
Bout=zeros(noj,1);                                                          % Placeholder for N x 1 amenity vector
Bout(Irsd)=EB;                                                              % Fill locations with positive residence employment. Those without receive theory-consistent zeros
HMC=zeros(noj,1);                                                           % Placeholder for N x 1 predicted workplace employment vector
HMC(Iwpl)=EHMC;                                                             % Fill locations with positive workplace employment. Those without receive theory-consistent zeros
HRC=zeros(noj,1);                                                           % Placeholder for N x 1 predicted residence employment vector
HRC(Irsd)=EHRC;                                                             % Fill locations with positive residence employment. Those without receive theory-consistent zeros
wout=zeros(noj,1);                                                          % Placeholder for N x 1 wage vector
wout(Iwpl)=Ewage;                                                           % Fill locations with positive workplace employment. Those without receive theory-consistent zeros

% We compute total worker income according 
% Amended code by GA to save memory
clear Epp_i Epp_j 
Epp_iji=Ephi_ij./((Ephi_i*ones(1,nwpl)));                                   % To compute the conditional commuting probability, we normalize bilateral commuting probabilitiy by the residential choice probability. The multiplication by the 1 x nwpl row vector of ones serves to replicate the residential choice probabilities nwpl times so that we obtain an nrsd x nwpl matrix that can be used to normalize the bilateral commuting probability matrix                   
EEWI = Epp_iji * Ewage  ;                                                   % We compute expected working income according to Eq. S.20  
ETWI = EEWI .* EHRT ;                                                       % We multiply expected worker income by residence empoyment to obtain total worker income
%%% Original code is outcommented

% Original outcommented code below
%{
Ephi_ji=Ephi_ij';                                                           % We transpose phi_ij to obtain a matrix that has workplaces in rows and residences in columns                                                                                         
Epp_iji=Ephi_ij./((Ephi_i*ones(1,nwpl)));                                   % To compute the conditional commuting probability, we normalize bilateral commuting probabilitiy by the residential choice probability. The multiplication by the 1 x nwpl row vector of ones serves to replicate the residential choice probabilities nwpl times so that we obtain an nrsd x nwpl matrix that can be used to normalize the bilateral commuting probability matrix                   
Epp_ijj=Ephi_ji./((Ephi_j*ones(1,nrsd)));                                   % We use the transposed matrix to compute the probability of living in i conditional on working in j 
temp1=Ewage.*EHMT;                                                          % We compute total worker income at workplace
temp2=Epp_ijj';                                                             % We transpose the transposed matrix so that we execute the matrix mutliplication 
temp3=(temp2*temp1);                                                        % We exectute the matrix multiplication of the conditional commuting probability matrix and the total worker income
vv=zeros(noj,1);                                                         % We create a n x 1 vector for the final outcome
vv(Irsd)=temp3;                                                          % We write the total worker income for blocks with positive residence employment into the n x1 vector, those with resident employment obtain a theory-consistent zero total worker income
%}
vv=zeros(noj,1);                                                         % We create a n x 1 vector for the final outcome
vv(Irsd)=ETWI;                                                           % We write the total worker income for blocks with positive residence employment into the n x1 vector, those with resident employment obtain a theory-consistent zero total worker income


% % Total worker income test;
% commute=Epp_ij.*HH;
% EWW=repmat(Ewage',nrsd,1);
% wcommute=commute.*EWW;
% test=zeros(noj,1);
% test(Irsd)=sum(wcommute,2);

% Commuting market access;
ECMA=d_ij_eps*(Ewage.^epsilon);                                             % We perform a matrix multiplication of the spatial weights matrix of dimension nfrom x nto and the nto x 1 vector of adjusted wages. This gives us the nfrom x 1 vector of (residential) commuting market access defined in Section S.3.1.2. 
CMA=zeros(noj,1);                                                           % Placeholder for full N x 1 vector
CMA(Irsd)=ECMA;                                                             % We fill locations with positive residence employment. The other locations remain with zero values

% Unconditional commuting probabilities;
ucprob=zeros(noj,noj);                                                      % Placeholder for unconditional commuting probabilities
ucprob(Irsd,Iwpl)=Epp_ij;                                                   % We fill routes with positive residence and workplace employment. The other routes remain with theory-consistent zeros


% Check equilibrium conditions;
% 
% test_i=sum(Epp_i);
% test_i
% test_j=EHMC-(Epp_iji'*EHRC);
% mean(test_j)
% test_i=EHRC-(Epp_ijj'*EHMC);
% mean(test_i)
% test_i=EHRC-(Epp_i.*HH);
% mean(test_i)
% test_j=EHMC-(Epp_j.*HH);
% mean(test_j)




% display('>>> Check Equilibrium Conditions <<<<');
% test_j=EHMT-(sum(Epp_ij)'*HH);
% mean(test_j)
% test_i=EHRT-(sum(Epp_ij')'.*HH);
% mean(test_i)


