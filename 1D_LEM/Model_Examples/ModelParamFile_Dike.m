%% Script Description
% Name: ModelParamFile_Dike
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Structured model variables (f,v,en,er,P,L) to run
% the 1D landscape evolution model from O'Hara et al. (2019), which allows 
% for different uplift and lithology parameters.
%
% *Script gives an example of an uplifting dike, defined as a vertical
%   layer with typical thicknesses and widths, and erodibility higher than 
%   the background. This shows the second of two ways in which vertical
%   beds can modeled.
%
% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.
    
%% File Parameters
f.filepath = pwd;
f.filename = 'LEM_Dike';

%% Model Parameters:
xmax = 1e4;
v.dx = 10;              
v.nx = xmax/v.dx;         
v.xRange = [0 xmax];
v.savePlot = 0;
v.allPlot = 2;
v.sinkHeight = 0;
v.landslide = 0;
v.critA = 30;
v.xr0 = (v.dx*v.nx)/2; % Initial ridge location (assume steady state).
v.minTime = 50e6;
v.ssRHeightPer = 1;
v.ssRPosPer = 1;
v.initRHeight = 1;
v.initRPos = 1;

%% End model Parameters:
en.st = 4e4;
en.tf = 20e6;
en.dt = 4e4;

%% Erosion model Parameters:
er.m = .5;
er.n = 1; % Model setup requires n to be 1.
er.h = 1.666;
er.ka = .5708;
er.u0 = 5e-5;
er.k = 4e-6;
er.D = .1;
er.Ps = 1;
er.Px = 1;

%% Perturbation Parameters:
P.runPerturbationModel = 0;
P.xP = 1.5217e+04;
P.R = 14000;
P.umax = .95;
P.b = 2;
P.uT = 1000;
P.dt = 10;
P.st = 10;
P.maxuT = 1000;

%% Lithology model Parameters:
L.runLithologyModel = 1;
L.plotLithology = 1;
L.kL = er.k/10;
L.DL = .1;
L.wL = 10000000;
L.xL = 2000;
L.thetaL = 90;
L.dL = 300;
L.Thickness = 100;

%% Uplift gradient parameters
U.runUpliftGradientModel = 0;
U.leftU = er.u0;
U.rightU = er.u0*10;

%% Run Model
addpath(genpath('./../'))
CreateAndRunModel(f,v,en,er,P,L,U)