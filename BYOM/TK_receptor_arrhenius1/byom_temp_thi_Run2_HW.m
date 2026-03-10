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
0.000	0.000	0.000	0.000	0.003	0.004	0.004
0.105	0.107	0.145	0.113	0.145	0.111	0.134
0.210	0.241	0.206	0.230	0.197	0.289	0.195
0.378	0.263	0.232	0.282	0.220	0.222	0.200
1.083	0.212	0.257	0.251	0.194	0.216	0.246
3.017	0.222	0.182	0.207	0.166	0.169	0.176
8.053	0.207	0.213	0.187	0.203	0.200	0.208

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
0.000	0.049
0.098	0.049
0.202	0.037
0.374	0.003
1.077	0.000
3.012	0.000
8.041	0.000


    ];

make_scen(4,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  


Temp1 = [ 0	2
0.000	26.1
0.042	26.0
0.083	25.7
0.125	25.1
0.167	23.9
0.208	23.1
0.250	22.9
0.292	22.6
0.333	22.2
0.375	22.2
0.417	22.3
0.458	22.6
0.500	22.7
0.542	22.7
0.583	22.6
0.625	22.6
0.667	22.6
0.708	22.8
0.750	23.0
0.792	23.3
0.833	23.4
0.875	24.3
0.917	25.7
0.958	26.3
1.000	26.3
1.042	26.3
1.083	25.7
1.125	25.3
1.167	24.9
1.208	24.5
1.250	24.0
1.292	23.3
1.333	22.7
1.375	22.2
1.417	22.3
1.458	22.5
1.500	22.6
1.542	22.5
1.583	22.7
1.625	23.1
1.667	23.6
1.708	24.0
1.750	24.3
1.792	25.0
1.833	26.0
1.875	26.6
1.917	27.0
1.958	27.1
2.000	27.1
2.042	27.0
2.083	26.3
2.125	25.7
2.167	25.0
2.208	24.0
2.250	23.5
2.292	23.5
2.333	23.5
2.375	23.6
2.417	23.5
2.458	23.0
2.500	22.8
2.542	22.8
2.583	22.5
2.625	22.4
2.667	22.3
2.708	22.5
2.750	22.8
2.792	23.6
2.833	24.3
2.875	25.2
2.917	25.9
2.958	25.9
3.000	26.0
3.042	26.3
3.083	26.7
3.125	26.2
3.167	25.4
3.208	24.4
3.250	23.4
3.292	22.5
3.333	21.7
3.375	21.0
3.417	20.2
3.458	19.6
3.500	19.3
3.542	19.1
3.583	19.0
3.625	19.0
3.667	19.1
3.708	19.7
3.750	20.7
3.792	23.1
3.833	24.3
3.875	25.3
3.917	25.9
3.958	26.4
4.000	26.9
4.042	27.0
4.083	26.5
4.125	24.8
4.167	24.2
4.208	24.0
4.250	23.7
4.292	23.3
4.333	22.6
4.375	22.1
4.417	21.6
4.458	21.4
4.500	21.3
4.542	21.4
4.583	21.5
4.625	21.8
4.667	22.1
4.708	22.4
4.750	22.9
4.792	23.6
4.833	23.8
4.875	25.0
4.917	25.8
4.958	26.2
5.000	26.1
5.042	25.3
5.083	24.0
5.125	23.0
5.167	22.2
5.208	21.7
5.250	21.2
5.292	20.7
5.333	20.3
5.375	19.8
5.417	19.5
5.458	19.2
5.500	18.9
5.542	18.7
5.583	18.5
5.625	18.5
5.667	18.4
5.708	18.9
5.750	19.7
5.792	20.5
5.833	21.1
5.875	22.7
5.917	24.1
5.958	24.8
6.000	24.9
6.042	24.3
6.083	23.8
6.125	23.2
6.167	22.6
6.208	20.9
6.250	18.4
6.292	16.2
6.333	14.5
6.375	13.3
6.417	12.4
6.458	11.7
6.500	11.1
6.542	10.8
6.583	10.5
6.625	10.3
6.667	10.2
6.708	10.7
6.750	11.7
6.792	13.1
6.833	14.3
6.875	15.2
6.917	15.7
6.958	16.1
7.000	16.4
7.042	16.3
7.083	15.9
7.125	15.5
7.167	14.8
7.208	14.2
7.250	13.8
7.292	13.5
7.333	13.3
7.375	13.3
7.417	13.4
7.458	13.5
7.500	13.6
7.542	13.8
7.583	13.8
7.625	13.9
7.667	13.9
7.708	14.1
7.750	14.7
7.792	15.7
7.833	17.1
7.875	18.8
7.917	19.8
7.958	21.0
8.000	19.4
8.042	18.3
8.083	16.7
8.125	14.8
8.167	13.6
8.208	13.2
8.250	12.9
8.292	12.8
8.333	12.7
8.375	12.6
8.417	12.6
8.458	12.7
8.500	12.7
8.542	12.6
8.583	12.6
8.625	12.5
8.667	12.4
8.708	12.4
8.750	12.8
8.792	13.9
8.833	15.1
8.875	15.8
8.917	16.0
8.958	15.2
9.000	15.8
9.042	15.2
9.083	14.5
9.125	13.8
9.167	13.2
9.208	13.0
9.250	12.8
9.292	12.7
9.333	12.5
9.375	12.4
9.417	12.4
9.458	12.2
9.500	12.2
9.542	12.2
9.583	12.2
9.625	12.3
9.667	12.5
9.708	12.7
9.750	13.0
9.792	13.7
9.833	15.1
9.875	16.2
9.917	16.6
9.958	17.0
10.000	17.0

];
      
make_scen(4, Temp1);  % linear extrapolation; Temp1 has header [0 2] => scenario 2
glo.T_scen = 2;        % tell derivatives/read_scen which scenario to use for temperature


% Create a table with nicer labels for the legends
Scenario = [1]; 
Label = {'Run1 Heatwave'};
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