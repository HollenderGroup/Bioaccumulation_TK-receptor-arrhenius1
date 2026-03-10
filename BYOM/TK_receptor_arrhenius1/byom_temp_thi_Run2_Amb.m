%% BYOM, byom_calanus_2016_onecomp.m 
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics.  
%
% *This script:* One-compartment TK of C2-naphthalene in _Calanus finmarchicus_.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 
%
%
% September 2022 
% Modifications by AMD to include (irrevirsible) receptor binding of an 
% antagonist. Here with the example of the antagonist thiacloprid (THI) 
% binding to the nicotinic-acetylcholine receptor (nAChR) as observed 
% in Gammarus pulex. in Raths et al. 202X


%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
% set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 1; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% Internal concentrations of THI in Gammarus pulex in [µmol/kg] 
DATA{3} = [0.5	1	1	1	1	1	1
0.000	0.000	0.000	0.000	0.003	0.004	0.003
0.103	0.097	0.101	0.108	0.061	0.051	0.076
0.209	0.140	0.098	0.140	0.110	0.098	0.144
0.377	0.143	0.207	0.200	0.186	0.207	0.120
1.082	0.126	0.134	0.076	0.161	0.174	0.157
3.016	0.142	0.197	0.174	0.129	0.156	0.154
8.044	NaN	NaN	NaN	0.147	0.152	0.155

];

% If needed, weights for individual measurements can be defined
% For this, uncommend the following line and specify your weights

% W{1} = 21 * ones(size(DATA{1})-1); % each point is a pooled sample of 21 animals

% In this data set, exposure was time-varying and reported as a series of
% concentrations over time. Here, the scenario is used as a linear forcing
% series (which has an analytical solution, and is thus much faster than
% the ODE version). Double time entries are used, which is more efficient,
% and probably more accurate.
% [µmol/L]
Cw1 = [ 0	1
0.000	0.036
0.099	0.035
0.203	0.034
0.376	0.003
1.078	0.000
3.014	0.000
8.042	0.000

    ];

make_scen(4,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  


Temp1 = [ 0	2
0.000	17.7
0.042	17.4
0.083	16.9
0.125	16.5
0.167	15.8
0.208	15.4
0.250	15.2
0.292	15.1
0.333	15.0
0.375	14.9
0.417	14.7
0.458	14.7
0.500	14.6
0.542	14.6
0.583	14.7
0.625	14.7
0.667	14.7
0.708	14.9
0.750	15.2
0.792	15.5
0.833	15.7
0.875	16.4
0.917	17.5
0.958	18.1
1.000	18.0
1.042	18.0
1.083	17.0
1.125	16.6
1.167	16.2
1.208	15.9
1.250	15.4
1.292	14.7
1.333	14.3
1.375	13.9
1.417	14.1
1.458	14.4
1.500	14.6
1.542	14.5
1.583	14.7
1.625	14.9
1.667	15.1
1.708	15.4
1.750	15.7
1.792	16.4
1.833	17.4
1.875	18.1
1.917	18.6
1.958	18.9
2.000	19.0
2.042	18.8
2.083	18.1
2.125	17.4
2.167	16.4
2.208	15.7
2.250	15.3
2.292	15.5
2.333	15.5
2.375	15.5
2.417	15.4
2.458	15.2
2.500	15.2
2.542	15.1
2.583	15.2
2.625	15.3
2.667	15.3
2.708	15.4
2.750	15.8
2.792	16.7
2.833	17.4
2.875	18.2
2.917	18.8
2.958	19.1
3.000	19.6
3.042	19.3
3.083	18.6
3.125	17.5
3.167	16.6
3.208	15.8
3.250	15.2
3.292	14.8
3.333	14.4
3.375	14.1
3.417	13.8
3.458	13.6
3.500	13.7
3.542	14.0
3.583	14.2
3.625	14.5
3.667	14.8
3.708	15.0
3.750	15.3
3.792	15.8
3.833	16.5
3.875	17.3
3.917	17.6
3.958	18.1
4.000	18.5
4.042	18.4
4.083	18.1
4.125	16.8
4.167	16.1
4.208	15.7
4.250	15.4
4.292	15.0
4.333	14.5
4.375	14.0
4.417	13.6
4.458	13.5
4.500	13.6
4.542	13.8
4.583	13.9
4.625	14.1
4.667	14.3
4.708	14.5
4.750	14.8
4.792	15.6
4.833	16.1
4.875	17.1
4.917	18.1
4.958	18.9
5.000	19.2
5.042	18.2
5.083	16.6
5.125	15.3
5.167	14.5
5.208	14.0
5.250	13.5
5.292	13.1
5.333	12.7
5.375	12.3
5.417	12.1
5.458	11.9
5.500	11.6
5.542	11.4
5.583	11.3
5.625	11.3
5.667	11.3
5.708	11.9
5.750	12.7
5.792	13.5
5.833	14.0
5.875	15.3
5.917	16.9
5.958	17.6
6.000	18.0
6.042	17.4
6.083	16.3
6.125	15.2
6.167	14.4
6.208	13.7
6.250	13.1
6.292	12.5
6.333	12.2
6.375	11.8
6.417	11.5
6.458	11.3
6.500	11.0
6.542	10.9
6.583	10.8
6.625	10.7
6.667	10.6
6.708	11.1
6.750	12.1
6.792	13.4
6.833	14.4
6.875	15.2
6.917	15.6
6.958	16.0
7.000	16.3
7.042	16.1
7.083	15.7
7.125	15.3
7.167	14.7
7.208	14.2
7.250	13.8
7.292	13.6
7.333	13.4
7.375	13.5
7.417	13.6
7.458	13.7
7.500	13.8
7.542	13.9
7.583	14.0
7.625	14.0
7.667	14.1
7.708	14.3
7.750	14.9
7.792	15.8
7.833	17.0
7.875	18.2
7.917	19.4
7.958	20.8
8.000	19.2
8.042	18.0
8.083	16.4
8.125	14.6
8.167	13.6
8.208	13.3
8.250	13.0
8.292	12.9
8.333	12.8
8.375	12.7
8.417	12.8
8.458	12.9
8.500	12.8
8.542	12.7
8.583	12.8
8.625	12.7
8.667	12.6
8.708	12.6
8.750	12.9
8.792	13.9
8.833	14.9
8.875	15.5
8.917	15.7
8.958	15.1
9.000	16.0
9.042	15.2
9.083	14.5
9.125	13.8
9.167	13.2
9.208	13.0
9.250	12.9
9.292	12.7
9.333	12.6
9.375	12.5
9.417	12.5
9.458	12.3
9.500	12.4
9.542	12.4
9.583	12.4
9.625	12.4
9.667	12.7
9.708	12.9
9.750	13.2
9.792	13.8
9.833	15.1
9.875	16.0
9.917	16.4
9.958	16.8
10.000	16.9

];
      
make_scen(4, Temp1);  % linear extrapolation; Temp1 has header [0 2] => scenario 2
glo.T_scen = 2;        % tell derivatives/read_scen which scenario to use for temperature


% Create a table with nicer labels for the legends
Scenario = [1]; 
Label = {'Run1 Ambient'};
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = Scenario; % scenarios (concentrations or identifiers)
X0mat(2,:) = 0;      % initial values state 1 (structure internal concentrations)
X0mat(3,:) = 0;      % initial values state 2 (receptor-antagonist complex concentration)
X0mat(4,:) = 0;      % initial values state 3 (total internal concentrations)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

glo.FMS = 0.01 ; % normalize for membrane protein content assuming a density of 1 (VMP/VS)

% syntax: par.name = [startvalue fit(0/1) minval maxval scale];
par.ke    = [3.2242         0       1      20       1];  % elimination rate constant, d-1
par.ku    = [8.5652         0       1      50       1];  % uptake rate constant, L/kg/d

par.TA_ke = [8336       0         1e3         5e6     1];   % Kelvin (Ea/R)
par.TA_ku = [7880       0         1e3         5e6     1];   % Kelvin (Ea/R)

par.kon   = [196.2      0   10     200      0];  % association of ligand-receptor complex
% % % % % %par.koff  = [0.00001  1     0 10 1];  % dissociation of ligand-receptor complex
par.B_MAX = [25      0   1      50       1];  % maximal binding capacity, µmol/kg


glo.Tref_C = 16;
glo.Tref_K = glo.Tref_C + 273.15;
glo.R_mod = 0;   % receptor model version (0 = none / not used) %placeholder for adding multiple model versions


%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = ['Concentration in structure compartment (',char(181),'mol/kg)'];
glo.ylab{2} = ['Concentration in membrane protein compartment (',char(181),'mol/kg)'];
glo.ylab{3} = ['Total internal concentration (',char(181),'mol/kg)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

glo.t = linspace(0, 10, 2000);   % dense
prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.
% par_out = calc_optim(par,opt_optim); % start the optimisation
% calc_and_plot(par_out,opt_plot); % call the plotting routine again to plot fits with CIs

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

opt_optim.type = 4; % optimisation method 1) simplex, 4) parameter-space explorer
opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 0; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.bw    = 0; % plot in black and white
opt_plot.cn    = 0; % if set to 1, connect model line to points (only for bw=1)
opt_plot.annot = 1; % annotations in sub-plot: text box with parameter estimates or overall legend
glo.useode     = 1; % use the analytical solution in simplefun.m (0) or the ODE solution in derivatives (1)
glo.stiff      = 0; % use ode45 with very strict tolerances

opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
% no plotting here; we'll immediately plot with CIs below

%% Plot results with confidence intervals
% The following code can be used to make a standard plot (the same as for
% the fits), but with confidence intervals. Options for confidence bounds
% on model curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 0; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time




%% Plot results (with CIs if available)

out_conf = [];              % default: no confidence intervals

opt_conf.type    = 3;       % keep if you WANT CIs from parspace explorer
opt_conf.lim_set = 0;
opt_conf.sens    = 0;

try
    out_conf = calc_conf(par_out,opt_conf);
catch ME
    % If no CI sample exists (common), just continue without CIs
    out_conf = [];
end

calc_and_plot(par_out,opt_plot,out_conf);  % <-- THIS makes the final plot