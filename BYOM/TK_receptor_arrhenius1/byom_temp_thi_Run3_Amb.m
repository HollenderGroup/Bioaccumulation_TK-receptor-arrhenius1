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
0.000	0.000	0.000	0.000	0.004	0.006	0.004
0.114	0.051	0.058	0.104	0.068	0.070	0.079
0.201	0.072	0.105	0.070	0.060	0.092	0.111
0.281	0.095	0.089	0.073	0.070	0.090	0.116
0.889	0.076	0.124	0.108	0.106	0.115	0.117
1.839	0.137	0.092	0.091	0.074	0.088	0.094
3.907	0.063	0.081	0.106	0.107	0.096	0.104

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
0.000	0.037
0.110	0.037
0.199	0.008
0.279	0.002
0.822	0.000
1.881	0.000
3.836	0.000
3.906	0.000

    ];

make_scen(4,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  


Temp1 = [ 0	2
0.000	17.1
0.042	16.8
0.083	16.4
0.125	16.2
0.167	16.0
0.208	15.8
0.250	15.7
0.292	15.6
0.333	15.5
0.375	15.4
0.417	15.2
0.458	14.9
0.500	14.7
0.542	14.4
0.583	14.1
0.625	14.1
0.667	14.2
0.708	14.6
0.750	15.1
0.792	15.8
0.833	17.0
0.875	18.2
0.917	18.3
0.958	18.3
1.000	18.1
1.042	17.1
1.083	16.1
1.125	15.3
1.167	14.7
1.208	14.2
1.250	13.6
1.292	13.6
1.333	13.7
1.375	13.9
1.417	14.1
1.458	14.0
1.500	14.0
1.542	13.9
1.583	13.9
1.625	13.8
1.667	13.8
1.708	14.1
1.750	14.4
1.792	14.8
1.833	15.2
1.875	15.6
1.917	16.1
1.958	17.2
2.000	17.5
2.042	16.9
2.083	16.2
2.125	15.5
2.167	14.7
2.208	14.1
2.250	13.5
2.292	13.1
2.333	13.1
2.375	13.4
2.417	13.7
2.458	13.8
2.500	13.9
2.542	14.0
2.583	14.0
2.625	14.0
2.667	14.1
2.708	14.3
2.750	14.7
2.792	15.2
2.833	15.8
2.875	16.1
2.917	16.6
2.958	17.2
3.000	17.0
3.042	16.5
3.083	15.9
3.125	14.9
3.167	14.3
3.208	14.2
3.250	14.1
3.292	14.0
3.333	14.0
3.375	14.0
3.417	13.8
3.458	13.5
3.500	13.5
3.542	13.5
3.583	13.4
3.625	13.4
3.667	13.5
3.708	13.8
3.750	14.4
3.792	15.4
3.833	15.9
3.875	15.6

];
      
make_scen(4, Temp1);  % linear extrapolation; Temp1 has header [0 2] => scenario 2
glo.T_scen = 2;        % tell derivatives/read_scen which scenario to use for temperature


% Create a table with nicer labels for the legends
Scenario = [1]; 
Label = {'50 ug/L'};
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

glo.t = linspace(0, 5, 2000);   % dense
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