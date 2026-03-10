%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files for Calanus2016, fit body residues with
% two-compartment TK. As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

C_S   = X(1); % state 1 (internal concentrations) at previous time point
C_MP  = X(2); % state 2 (receptor-agonist complex concentration) at previous time point

% these concentrations are both expressed on volume basis.

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

FMS = glo.FMS;  % normalize for membrane protein content assuming a density of 1 (VMP/VS)



ke_ref = par.ke(1);     % elimination rate constant at reference temperature, d-1
ku_ref = par.ku(1);     % uptake rate constant at reference temperature, L/kg/d


TA_ke  = par.TA_ke(1); %Arrhenius temperature ke, K
TA_ku  = par.TA_ku(1); %Arrhenius temperature ku, K


kon   = par.kon(1);   % association of ligand-receptor complex
% %koff  = par.koff(1);  % dissociation of ligand-receptor complex
B_MAX = par.B_MAX(1); % maximal binding capacity, Âµmol/kg 


%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed).

% ---- Exposure concentration forcing (°C) ----
if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)
end

% ---- Temperature forcing (°C) ----
T_C = read_scen(-1, glo.T_scen, t, glo);   % scenario index stored in glo.T_scen
T_K = T_C + 273.15;
Tref = glo.Tref_K;

% T_K  = glo.Tref_K;   % force constant temperature at reference
% Tref = glo.Tref_K;

%% Calculate the derivatives
% This is the actual model, specified as a system of two ODEs:

% Calculating the change of ligand concentration in membrane protein
% fraction with either (1) the Michaelis-Menten kinetics or (2) the second
% order kinetics. 

% calculate temperature dependent rate constants
ku = ku_ref * exp( TA_ku * (1/Tref - 1/T_K) );
ke = ke_ref * exp( TA_ke * (1/Tref - 1/T_K) );

% Calculating the change of the ligand concentration in the membrane
% protein compartment.
dC_MP =  kon * C_S * (B_MAX - C_MP) ; % second order kinetics

% Calculating the change of the ligand concentration in the structure
% compartment.
dC_S = ku * c - ke * C_S - dC_MP * FMS; % first order bioconcentration

%Calculating the change of the total ligand concentration in the organism.
dC_tot = dC_S  + dC_MP * FMS ;

dX = [dC_S; dC_MP; dC_tot]; % collect all derivatives in one vector dX, ***********dC_MP;
