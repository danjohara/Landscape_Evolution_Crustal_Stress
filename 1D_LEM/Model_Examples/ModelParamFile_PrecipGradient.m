%% Script Description
% Name: ExampleModelParamFile_UpliftPerturbation
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Structured model variables (f,v,en,er,P,L) to run
% the 1D landscape evolution model from O'Hara et al. (2019), which allows 
% for different uplift and lithology parameters.
%
% *Script gives an example of an uplift pertrubation similar to O'Hara 
%   et al (2019), subject to a precipitation gradient.
%
% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.
    
%% File Parameters
f.filepath = pwd;
f.filename = 'LEM_PGrad';

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
v.minTime = 40e6;
v.ssRHeightPer = 1;
v.ssRPosPer = 1;
v.initRHeight = 1;
v.initRPos = 1;

%% End model Parameters:
en.st = 4e4;
en.tf = 30e6;
en.dt = 4e4;

%% Erosion model Parameters:
er.m = .5;
er.n = 1; % Model setup requires n to be 1.
er.h = 1.666;
er.ka = .5708;
er.u0 = 5e-5;
er.k = 4e-6;
er.D = .1;
er.Px = v.xRange(1):v.dx:v.xRange(2);

highPrecLoc = 4000;
er.Ps = exp(-.5*((er.Px-highPrecLoc)/4000).^2);
t2 = exp(-.5*((er.Px-highPrecLoc)/1000).^2);
er.Ps(er.Px>highPrecLoc) = t2(er.Px>highPrecLoc);
er.Ps(er.Ps<.5) = .5;

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
L.runLithologyModel = 0;
L.plotLithology = 1;
L.kL = 5e-6;
L.DL = .01;
L.wL = 30000;
L.xL = 30000;
L.thetaL = 10;
L.dL = 4000;
L.Thickness = 500;

%% Uplift gradient parameters
U.runUpliftGradientModel = 0;
U.leftU = er.u0;
U.rightU = er.u0;

%% Run Model
addpath(genpath('./../'))
CreateAndRunModel(f,v,en,er,P,L,U)