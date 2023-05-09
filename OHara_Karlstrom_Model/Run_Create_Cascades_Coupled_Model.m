%% 
% Name: Run_Create_Cascades_Coupled_Model
% Author: Daniel O'Hara
% Date: 05/09/2023
%
% Description: Script to run the coupled landscape evolution - crustal
% stress model for a Cascades-like topography and precipitation, creating
% the figures shown in O'Hara & Karlstrom (in review).
%
% Reference:
% O'Hara, D., and Karlstrom, L., in review, Distributed volcano erosion 
%   patterns imply coupled volcanism and regional climate in the Cascades 
%   arc, U.S.A. Frontiers in Earth Science.

%% Add path
clear
addpath(genpath('.\..\'))

%% Run landscape evolution model
Run_Topography_Model

%% Create crustal stress plots
Plot_Model_Stresses

%% Make movie of the model
Make_Model_Movie