%% Script Description
% Name: LEM_Stresses_Dike
% Author: Daniel O'Hara
% Date: 05/09/2023
%
% Description: Script to run the coupled landscape evolution - crustal 
% stress model from O'Hara & Karlstrom (in review).
%
% *Script gives an example of an uplifting dike, defined as a vertical
%   layer with typical thicknesses and widths, and erodibility higher than 
%   the background. This shows the second of two ways in which vertical
%   beds can modeled.
%
% References:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.
%
% O'Hara, D., and Karlstrom, L., in review, Distributed volcano erosion 
%   patterns imply coupled volcanism and regional climate in the Cascades 
%   arc, U.S.A. Frontiers in Earth Science.

%% Run LEM First
addpath(genpath('./../../'))
cd Example_Results\
ModelParamFile_Dike
clear
close

%% Setup
modFile = 'LEM_Dike.mat';
saveFol = [pwd,'\'];
analyzeSteps = [1,26,101,501];
saveName = 'Dike_Step_%d';
scaleHeight = 2000;

inputs.G = .2e11;
inputs.lambda = 0;
inputs.rho = 2000;

inputs.dz = 10;
inputs.maxDepth = 4000;
inputs.minDepth = 10;

inputs.sigmaP_scale = 300;
inputs.sigmaP_cutByX = 50;
inputs.sigmaP_cutByZ = 100;
inputs.magScale1 = [-.01,.1];
inputs.magScale2 = [0,.5];

inputs.kernel_radius = 1e4;
inputs.Analysis = 1;

inputs.VectorWidth = 2;
inputs.makePlots = 1;
tmp = [2000:2000:8000]';
inputs.sourceLocs = [tmp,ones(size(tmp))*inputs.maxDepth-50];

%% Load Topgraphy
load(modFile)
inputs.x = v.x;
inputs.dx = abs(v.x(2)-v.x(1));

%% Determine Max Height
maxHeight = -Inf;
sF = 1;
for i = 1:length(analyzeSteps)
    tt = fv_all.zt(analyzeSteps(i),:);
    if i == 1
        sF = scaleHeight/max(tt(:));
    end
    tt = tt*sF;

    maxHeight = max([maxHeight,max(tt)]);
end

inputs.topoMaxHeight = maxHeight;

%% Run Analyses
scaleFactor = 1;
for i = 1:length(analyzeSteps)
    inputs.Topo = fv_all.zt(analyzeSteps(i),:);

    if i==1
        scaleFactor = scaleHeight/max(inputs.Topo(:));
    end
    inputs.Topo = inputs.Topo*scaleFactor;
    inputs.plotTitle = sprintf('T = %.2f Myr',fv_all.ts(analyzeSteps(i))/1e6);

    sigmas = Calculate_Crustal_Stress_2D(inputs);

    save([saveFol,sprintf(saveName,analyzeSteps(i)),'.mat'])

    orient landscape
    set(gcf,'PaperType','tabloid')
    set(gcf,'Renderer','opengl')
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_trajectories_bitmap.fig'],'fig');
    print([saveFol,sprintf(saveName,analyzeSteps(i)),'_trajectories_bitmap'],'-dpdf','-r300');
    close

    orient landscape
    set(gcf,'PaperType','tabloid')
    set(gcf,'Renderer','painters')
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_trajectories_vector.fig'],'fig');
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_trajectories_vector.pdf'],'pdf');
    close

    orient landscape
    set(gcf,'PaperType','tabloid')
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_trajectories.fig'],'fig');
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_trajectories.pdf'],'pdf');
    close

    orient landscape
    set(gcf,'PaperType','tabloid')
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_sigma1_topo.fig'],'fig');
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_sigma1_topo.pdf'],'pdf');
    close

    orient landscape
    set(gcf,'PaperType','tabloid')
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_sigmaX_sigmaZ.fig'],'fig');
    saveas(gcf,[saveFol,sprintf(saveName,analyzeSteps(i)),'_sigmaX_sigmaZ.pdf'],'pdf');
    close
end