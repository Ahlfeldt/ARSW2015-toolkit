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
%%% Then, it uses solved transformed wages to recoversadjusted productivity 

function [vvout] = expincome(obsdata,distvar,noj,wage,B);
    % This program uses the following inputs
        % obsdata refers to an n by 4 object that contains a vector of four variables,
            % floor space prices, worplace employment, residence employment, and area        
        % distvar is a n by n matrix of bilateral travel cost (km or time or euros)  
        % noj is a scalar that defines the number of observations n     
        % wage is a vector of guesses of ADJUSTED wages solved by comegaoptC.m
            % (notice, comegaoptO delivers TRANSFORMED wages)
        % B is a vector adjusted amenities solved by camen.m
    % This program produces the following outputs
        % vvout is total worker income (expected worker income times residence
        % employment
    % The names of the inputs need to correspond to objects that exist in
    % the workspace. The names of the outputs can be freely chosen

% Declaring scalars as globals so they can be accessed and modified by any function
global alpha beta kappa epsilon;

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
Ewage = wage(Iwpl);                                                         % Generates a vector of adjusted wages only containing observations with positive workplace employment
EB = B(Irsd);                                                               % Generates a vector of adjusted amenities only containing obserations with positive residence employment
EKM = K(Iwpl);                                                              % Generates a vector of geographic area only only containing observations with positive workplace employment                                                       
EKR = K(Irsd);                                                              % Generates a vector of geographic area only only containing observations with positive residence employment    
HH=sum(EHMT);                                                               % Compute total workplace employment

% Commuting costs;
Edistvar=distvar(Irsd,Iwpl);                                                % Generates a smaller travel time matrix that connects only observations with positive workplace employment to observations with residence emplyoment (it has nto by nfrom osbervations)
clear distvar                                                               % Notice that this time, we are having residence in rows and workplace in columsn, so we have an nfrom x nto matrix
d_ij_eps=exp(-epsilon.*kappa.*Edistvar);                                    % Generates a spatial weight matrix with weights that exponentially falls in bilateral travel time between i and j at rate epsilon*kappa
clear Edistvar                                                               % clear variable to save memory 

% Computing bilateral and conditional commuting probabilitieis in Eq. 4/5
% Workers live in i (rows = residences) and work in j (cols = workplaces);

%%% New GA code starts %%%
clear cprob06 % not used here can be dropped if in data
% Assuming the rest of the setup is correct and the variables have the expected sizes
% Preallocate matrix for Ephi_ij
Ephi_ij = zeros(nrsd, nwpl);
% Compute Ephi_ij element-wise
EQT=QT(Irsd);
for i = 1:nrsd
    for j = 1:nwpl
        % Assuming QT(Irsd(i)), EB(i), Ewage(j), and d_ij_eps(i, j) are scalar values
        EQQ_element = EQT(i)^(-(1-beta) * epsilon);
        EBB_element = EB(i)^epsilon;
        EWW_element = Ewage(j)^epsilon;
        Ephi_ij(i, j) = d_ij_eps(i, j) * EQQ_element * EBB_element * EWW_element;  % Correct as is
    end
end
% Clear variables to free up memory
clear d_ij_eps
%%% New GA code ends %%%
% Original outcommented code below
%{
EQQ=repmat(QT(Irsd),1,nwpl);                                                % We assign residence floor space prices to all bilaterals. We take the vector of floor space prices, generate a new vector only for observations with positive residence employment, and replicates it nwpl times to generate a nrsd x nwpl matrix
EQQ=EQQ.^(-(1-beta).*epsilon);                                              % Generate the first component of bilateral commuting probabilities in Eq. 4: floor space prices to the power of (1-beta) x -epsilon 
EBB=repmat(EB,1,nwpl);                                                      % We assign residence amenities to all bilaterals. Takes the vector of amenities for observations with positive residence emplyoment and replicates it to generate an nrds x nwpl matrix
EBB=EBB.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities: B to the power of epsilon 
EWW=repmat(Ewage',nrsd,1);                                                  % We assign workplace wages to all bilaterals. We take the nwpl x 1 vector of adjusted wages, transpose it into an 1 x nwpl vector and replicate in nrsd rows to obtain an nrsd x nwpl matrix
EWW=EWW.^epsilon;                                                           % Generates the next component of bilateral commuting probabilities, wage to the power of epsilon
Ephi_ij=d_ij_eps.*EBB.*EQQ.*EWW;                                            % Now we use the above components to compute the numerator of bilateral commuting probabilities in Eq. 4, phi_ij
clear EQQ EWW EBB d_ij_eps
%}

Ephi=sum(sum(Ephi_ij));                                                     % We sum over all bilaterals to obtain the denominator in Eq. 4. sum(Ephi_ij) retuns column sums, so we sum over the column sum to get the sum over all bilaterals
Ephi_i=sum(Ephi_ij'); Ephi_i=Ephi_i';                                       % We compute the numerators in Eq. 5: We transpose the phi_ij matrix so that residences are in columns. The row sum then returns sum(phi_ij) by i as a 1 x nrsd row vector. We transpose this vector to obtain a nrsd x 1 column vector.
Ephi_j=sum(Ephi_ij); Ephi_j=Ephi_j';                                        % We compute the numerators in Eq. 5: Summing over the columns of the nrsd x nwpl phi_ij matrix, returns sum(phi_ij) by j. We transpose the 1 x nwpl row vector to obtain a nwpl x 1 column vector.
Epp_ij=Ephi_ij./Ephi;                                                       % We divide the numerator by the denominator to obtain the commuting probabilities                                                  
Epp_i=Ephi_i./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a residence
Epp_j=Ephi_j./Ephi;                                                         % We divide the numerator in Eq. 5 by the denominator to obtain the probability that workers choose a workplace
EHRC=Epp_i.*HH;                                                             % Compute residence employment predicted commuting probabilities (should be equal to observed reidence employment)
EHMC=Epp_j.*HH;                                                             % Compute workplace employment predicted commuting probabilities (should be equal to observed reidence employment)

% We compute total worker income according 
% Amended code by GA to save memory
clear Ephi Epp_i Epp_j Epp_ij
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
vvout=zeros(noj,1);                                                         % We create a n x 1 vector for the final outcome
vvout(Irsd)=temp3;                                                          % We write the total worker income for blocks with positive residence employment into the n x1 vector, those with resident employment obtain a theory-consistent zero total worker income
%}
vvout=zeros(noj,1);                                                         % We create a n x 1 vector for the final outcome
vvout(Irsd)=ETWI;                                                           % We write the total worker income for blocks with positive residence employment into the n x1 vector, those with resident employment obtain a theory-consistent zero total worker income

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


