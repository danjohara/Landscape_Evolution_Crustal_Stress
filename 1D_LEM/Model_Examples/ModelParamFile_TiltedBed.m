%% Script Description
% Name: ExampleModelParamFile_TiltedBed
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
f.filename = 'Example_TiltedBed3';

%% Model Parameters:
xmax = 6e4;
v.dx = 50;              
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
en.st = 1e5;
en.tf = 15e6;
en.dt = 1e5;

%% Erosion model Parameters:
er.m = .5;
er.n = 1; % Model setup requires n to be 1.
er.h = 1.666;
er.ka = .5708;
er.u0 = 1e-4;
er.k = 1e-6;
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
L.kL = 1e-7;
L.DL = .01;
L.wL = 1000;
L.xL = 20000;
L.thetaL = 80;
L.dL = 500;
L.Thickness = 250;

%% Uplift gradient parameters
U.runUpliftGradientModel = 0;
U.leftU = er.u0;
U.rightU = er.u0*10;

%% Run Model
addpath(genpath('./../'))
CreateAndRunModel(f,v,en,er,P,L,U)