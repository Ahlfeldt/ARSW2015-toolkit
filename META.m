%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master MATLAB programme file for the toolkit for    %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: GMA, 01/2ß24                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program is not part of the orginal replication directory       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please the following toolboxes
    % Global optimization
    % Statistics and Machine learning 
    % Mapping

% Root folder of working directory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User 1 GA desktop
% User 2 GA notebook
% User 3 GA Server
% User 4 AH
user = 3                                                                    % Define users
if user==1;                                                                 % Root directory of primary user
    cd 'D:/Dropbox/_HUB_HerreraA/Course/Repository/Replication_Directories/ARSW2015-Teaching/matlab';             
elseif user==2;                                                             % Root directory of secondary user. Add more users if necessary
    cd 'C:/Users/gabri/Dropbox/_HUB_HerreraA/Course/Repository/Replication_Directories/ARSW2015-Teaching/matlab';
elseif user==3;                                                             % Root directory of secondary user. Add more users if necessary
    cd 'E:/Data/Gabriel/Ahfleldt/_QSE/ARSW/FinalBerlin/matlab';
elseif user==4;
    cd 'C:/Users/andre/Dropbox/_HUB_HerreraA/Course/Repository/Replication_Directories/ARSW2015-Teaching/matlab';
end;
addpath('section6/optimepsilon')                                            % Adding path to additional do files and functions that solve for transformed wages and estimate epsilon
addpath('section6/calibration')                                             % Adding path to additional do files and functions that solve for transformed wages and estimate epsilon
addpath('section6/exogcftual')                                              % Adding path to additional do files and functions that solve for transformed wages and estimate epsilon
addpath('section7/counterfactual')                                              % Adding path to additional do files and functions that solve for transformed wages and estimate epsilon


% Data preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepdata_TD86
    % generates preddata_TD86 data set for one-step estimation. No need to
    % exececute unless you want to add more variables from the replication
    % directory. Designed to minimize memory requirement for execution on
    % desktop computers
% predata_TD
    % generates preddata_TD containing 2006 data for quantification and
    % inversion. No need to exececute unless you want to add more variables
    % from the replication directory. Designed to minimize memory 
    % requirement for execution on desktop computers

% Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimepsilon_TD86 ;
    % solves for 1986 transformed wages and finds epsilon that results in
    % adjusted 1986 Bezirke wages with the same log variance as observed in
    % data. Generates maps of transformed and adjusted wages.
    % Uses 
        % section6/comegaopt.m: Finds the transformed wages for given
            % residence and workplace emplyoment
        % section6/cdensityoptren.m: Computes objective function that
            % search algorithm minimizes to find epsilon
        
% Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calcal_TD ;
    % Recovers structural fundamentals from observed endogenous variables
        % using the structure of the model
    % Uses a sequential procedure that is didactically useful
    % uses 
        % comegaoptC.m to solve for adjusted wages and recover adjusted
            % productivities
        % camen.m to recoved ajusted amenities
        % calcal_adj_TD.m to rescale amenities to rationalize city size
            % (new program)
        % expincome.m to recover total income by residence
        % cdensity.m to recover cdensity 
        
% Counterfactuals with exogenous fundamentals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cftualprep_TD ;
    % We prepare the data set for the counterfactuals by running a
        % simlutaneous inversion for amenities and productivities along with
        % solving for wages
        % These amenities and productivities ensure that the model recovers 
        % initial values of endogneous variables
    % We compare recovered productivities and amenities to those recovered
        % using the sequential procedure called by calcal_TD.m
cftualexog_TD  ;
    % We change primitives and solve for counterfactual outcomes under
        % exogenous fundamentals
        % In the closed-city model
    % We perform illustrative counterfactuals
        % Banning road-traffic
        % Improving productivity of East Berlin
    % Uses
        % smodexog.m to solve for endogenous variables for given primitives
        
% Counterfactuals with endogenous fundamentals %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cftualprep_end_TD ;    
    % We prepare the data set for the counterfactuals by running a
        % simlutaneous inversion for amenities and productivities along with
        % solving for wages
    % Unlike in in the case with exogenous fundamentals
        % we now use the GMM estimates of the full model in the
        % quantification of the model
    % We decompose the amenities and productivities in an ednogenous and an
        % exogenous part
    % We illustrate how after accounting for endogenous components,
        % fundamentals are "flat" (not much correlated with CBD distance)
    % We compare productivities and amenities under the GMM parametrization
        % to those obtained under parameter values from Section 6
cftualendog_HHbar_TD;
    % We conduct illustrative counterfactuals with endogenous agglomeration 
    % forces and constant total employment
        % To this end we hold population constant and let utility adjust 
        % These counterfactuals correspond to the closed-city case
  % We conduct two illustratve counterfactuals
        % A complete car ban: We solve for the counterfactual equilibrium
            % using public transit travel times
        % An upgrade to fundamental productivity in East Berlin: We
            % increase a with the boundaries of former East Berlin
   % For both counterfactuals, we provide a comparison to the case with
        % exogenous fundamentals        
cftualendog_Ubar_TD;
    % We conduct illustrative counterfactuals with endogenous agglomeration 
    % forces and constant reservation utility
        % To this end we adjust city population 
        % These counterfactuals correspond to the open-city case
    % We conduct two illustratve counterfactuals
        % A complete car ban: We solve for the counterfactual equilibrium
            % using public transit travel times
        % An upgrade to fundamental productivity in East Berlin: We
            % increase a with the boundaries of former East Berlin
    % For both counterfactuals, we provide a comparison to the 
        % open-city case