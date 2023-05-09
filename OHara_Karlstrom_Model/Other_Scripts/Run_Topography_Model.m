%% 
% Name: Run_Topography_Model
% Author: Daniel O'Hara
% Date: 05/09/2023
%
% Description: Script to create the landscape evolution model component for  
% the coupled landscape evolution - crustal stress model for a 
% Cascades-like topography and precipitation shown in O'Hara & Karlstrom 
% (in review).
%
% Reference:
% O'Hara, D., and Karlstrom, L., in review, Distributed volcano erosion 
%   patterns imply coupled volcanism and regional climate in the Cascades 
%   arc, U.S.A. Frontiers in Earth Science.

%% Setup
useModFol = '.\..\';
addpath(genpath(useModFol))

meanPFile = 'Cascades_meanZP.mat';

%% Base Model Setup
% Main Model Parameters:
xmax = 8e4;                 
v.dx = 20;                  % Grid resolution
v.nx = xmax/v.dx;           % Number of horizontal nodes.
v.xRange = [0 xmax];        % Horizontal domain size.
v.savePlot = 0;             % Flag to save landscape plots
v.allPlot = 2;              % Flag for which landscape components to plot.
v.sinkHeight = 0;           % Flag for how to treat sinks.
v.landslide = 0;            % Flag for landsliding (slope constraint).
v.critA = 30;               % Critical drainage area for channelization
v.xr0 = (v.dx*v.nx)/2;      % Initial ridge location (assume steady state).
v.minTime = 40e6;           % Minimum model runtime.
v.ssRHeightPer = 1;         % Percentile height for steady-state (ignored
v.ssRPosPer = 1;
v.initRHeight = 1;
v.initRPos = 1;

% Lithology model Parameters:
L.runLithologyModel = 0;    % Flag to run the lithology model.

% Uplift gradient parameters
U.runUpliftGradientModel = 0;   % Flag to run uplift gradient model.

% Base Erosion Parameters
er.m = .5;      % Drainage area expononent.
er.n = 1;       % Slope exponent, model setup requires n to be 1.
er.h = 1.666;   % Hack's Law exponenent.
er.ka = .5708;  % Hack's Law coefficient.
er.u0 = 1e-4;   % Tectonic uplift rate.
er.k = 6e-6;    % Rock erodibility.
er.D = 1e-1;    % Soil diffusivity.
er.Ps = 1;      % Precipitation scaling.

%% Perturbation Model
% Run uplift perturbation model for 2.6 Myr.
P.runPerturbationModel = 1;     % Flag for uplift perturbation model.
P.xP = 4.7e4;                   % Location of uplift locus.
P.R = 20e3;                     % Scale-length (radius) of uplift parabola.
P.umax = 7e-4;                  % Maximum perturbation uplift.
P.b = 2;                        % Exponent for parabola shape.
P.uT = 2.6e6;                   % Pertubation model runtime.
P.dt = 1e4;                     % Pertubation model timestep.
P.st = 1e4;                     % Model save steps.
P.maxuT = 2.6e6;                % Model time of uplift perturbation.

en.st = 1e4;                    % Post-uplift model savetime
en.tf = 3e4;                    % Post-uplift model runtime
en.dt = 1e4;                    % Post-uplift model timestep

%% Run Perturbation Model
f.filepath = '.\LEM_Model_Results\';
f.filename = 'L1_Uplift';

fv = CreateAndRunModel(f,v,en,er,P,L,U);

%% Create Climate-Perturbation Model
% String models together,
v.prevFv= fv;

% Update perturbation parameters.
P.runPerturbationModel = 1;
P.uT = 100e6;
P.maxuT = 100e6;
P.dt = 1e5;
P.st = 1e5;

% Reset time
en.st = 1e5;  
en.tf = 1e5; % Set to 100e6 if not running perturbation, 5e5 otherwise.
en.dt = 1e5;

% Update erosion params for precitation gradient
load(meanPFile)
tmpX = fv.v{2}.x;
tmpZ = fv.zt(end,:);
ii = find(tmpZ==max(tmpZ));
tmpX = tmpX-tmpX(ii);

modP = interp1(newPArcDists,meanMPR,tmpX);
modP = modP./max(modP);
er.Ps = modP;                               % Precipitation scaling
er.Px = v.xRange(1):v.dx:v.xRange(2);       % Precipitation x-values

%% Run Climate-Perturbation Model
f.filename = 'L2_Precip';
fv = CreateAndRunModel(f,v,en,er,P,L,U);