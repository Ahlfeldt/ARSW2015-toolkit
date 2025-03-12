%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB programme file for the toolkit for           %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: SJR, 02/2015                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This version has been commented                                     %%%
%%% The syntax has been updated to refer to Epsi (instead of theta) as  %%%
%%% as the Frechet shape parameter (consistent with notations in the    %%%
%%% paper                                                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = cdensityoptren(Epsi)    
    % Epsi is the seet of parameters patternsearch is attempting to
        % optimize. The initial value of Epsi is set to Epsi0 (set to 4 in optimespsilon.m). Patternsearch
        % iteratively updates the value of Epsi as it searches for the minimum
        % of the objective function. At each iteration, patternsearch calls
        % cdensityoptren with the current value of Epsi
    % f is the value of the objective function that is being retuned by the cdensityoptren function and then minimized by the patternsearch function in optimepsilon  

global epsilon kappaeps lambda delta eta rho alpha beta ftD eps Theta0;
global LD WD fwestd fwestr bzk86dw bzk86rw;
global obsvar36 tt36 nobs36 omega_36 omega36 wage36 a36 b36;
global obsvar06 tt06 nobs06 omega_06 omega06 wage06 a06 b06;
global obsvar86dw tt86dw nobs86dw omega_86dw omega86dw wage86dw a86dw b86dw;
global obsvar86rw tt86rw nobs86rw omega_86rw omega86rw wage86rw a86rw b86rw;
global IICBDdw IICBDrw IIBZKdw IIBZKrw;
global varlwdata varlwmod lbwdata lbwmod;

epsilon=Epsi(1);                                                            % We define epsilon as the value of the parameter (Epsi) over which the patternsearch algorithm is searching

empwpl86rw=obsvar86rw(:,2);                                                 % obsvar86rw includes floor space prices, workplace employment, residence employment, and area for West Berlin in 1986

% ***************;
% **** WAGES ****;
% ***************;

wage86rw=omega86rw.^(1./epsilon);                                           % Reverse transformation to get back from transformed wages to adjusted wages

wage86rw(wage86rw>0)=wage86rw(wage86rw>0)./geomean(wage86rw(wage86rw>0));   % Normalize wages by geometric mean within observations with positive wages

wbill86rw=wage86rw.*empwpl86rw;                                             % Generate 1986 wage bill as the product of wage and employment

bzknum86rw=grpstats(wbill86rw,bzk86rw,'numel');                             % Computes the number of elements (blocks) within each Bezirk (denoted by bryk86rw)
bzkwbill86rw=grpstats(wbill86rw,bzk86rw,'mean');                            % Computes the mean wage bill (blocks) within each Bezirk (denoted by bryk86rw)
bzkwbill86rw=bzkwbill86rw.*bzknum86rw;                                      % Multiply Bezirke mean wage bill by number of elements of Bezirk to get Bezirke wage bill
bzkempw86rw=grpstats(empwpl86rw,bzk86rw,'mean');                            % Computes the mean workplace employment (blocks) within each Bezirk (denoted by bryk86rw)
bzkempw86rw=bzkempw86rw.*bzknum86rw;                                        % Multiply Bezirke mean workplace employment by number of elements of Bezirk to get Bezirke wage bill
bzkwage86rw=bzkwbill86rw./bzkempw86rw;                                      % Divide total Bezirke wage bill by total Bezirke employment to get Bezirke wage

lbwmod=log(bzkwage86rw);                                                    % Compute log Bezirke wage
lbwmod=lbwmod-mean(lbwmod);                                                 % De-mean log-Bezirke wage
varlwmod=var(lbwmod);                                                       % We compute the variance of log wages

% ***************************;
% **** Moment Conditions ****;
% ***************************;

LD=size(varlwmod,1);                                                        % Recover the number of rows of the variance of log wages vector (here it is 1)                                                  
WD=eye(LD);                                                                 % Creates an identity matrix of the n dimension of the log wages variance vector (here just 1)
ftD=[varlwmod-varlwdata];                                                   % Comparing variance of log bezirke wages in model and data (here, stored as a scalar, could be vector)

ftt=ftD' * WD * ftD;                                                        % Here we do a matrix multiplication of the transposed vector of errors in log variances and the vector itself (the multiplication by the identity matrix WD is not consequentia). The result is the dot product which corresponds to the sum of squares; here, the result is the residial sum of squares.
f = ftt.*(10.^6);                                                           % We scale it up since the RSS are very small and we want to avoid numerical issues

% Notice that this moment condition corresponds to equation S.64 the
% supplement. This moment condition says
% E((1/epsilon)^2 ln(omega)^2-(sigma_(ln(w)))^2)=0. This notation makes
% explicit that we can use the moment condition to identify epsilon. In the
% code above we have created the adjusted wage wage86rw from 
% transformed wage, omega, in line 36. Therefore, epsilon no longer shows up
% when we compute the RSS between adjusted wages in the model and observed
% wages in data. Notice further that E((1/epsilon)^2 ln(omega)^2) is the
% variance of transformed wages since omega has a mean of one and, hence,
% ln(omega)=0.

