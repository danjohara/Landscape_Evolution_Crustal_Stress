%% 
% Name: Plot_Model_Stresses
% Author: Daniel O'Hara
% Date: 05/09/2023
%
% Description: Script to plot the crustal stress componenets for the 
% coupled landscape evolution - crustal stress model for a Cascades-like 
% topography and precipitation shown in O'Hara & Karlstrom (in review).
%
% Reference:
% O'Hara, D., and Karlstrom, L., in review, Distributed volcano erosion 
%   patterns imply coupled volcanism and regional climate in the Cascades 
%   arc, U.S.A. Frontiers in Earth Science.

%% Setup
modFile = '.\LEM_Model_Results\L2_Precip.mat';
scriptPath = '.\..\';
saveFol = '.\Crustal_Stress_Plots\';
meanPFile = '.\Cascades_Data\Cascades_meanZP.mat';

stepTitles = {'Initial','Uplift','Precip'};

addpath(genpath(scriptPath))

%% Load File and Collect Model Values
load(modFile)

v = v_uplift;
inputs.dz = v.dx*5; % Vertical grid resolution
inputs.x = v.x;     % Horizontal coordinates.
inputs.dx = abs(v.x(2)-v.x(1)); % Horizontal grid resolution.

useZs = [v.prevFv.zt(1,:);v.prevFv.zt(end,:);fv_all.zt(end,:)];
useTs = [v.prevFv.ts(1);v.prevFv.ts(end);fv_all.ts(end)+v.prevFv.ts(end)];

%% Other Parameter values
inputs.G = 10e9;
inputs.lambda = 20e9;
inputs.rho = 2700;

inputs.maxDepth = 10050;
inputs.minDepth = 10;

inputs.sigmaP_scale = 2000;
inputs.sigmaP_cutByX = 400;
inputs.sigmaP_cutByZ = 30;
inputs.magScale1 = [-.01,.1];
inputs.magScale2 = [0,.5];

inputs.kernel_radius = 2e4;
inputs.Analysis = 1;

inputs.VectorWidth = 2;
inputs.makePlots = 1;

tmp = [30e3:5e3:50e3]';
inputs.sourceLocs = [tmp,ones(size(tmp))*inputs.maxDepth-50];

%% Make forcing figures
vUplift = v.prevFv.v{1};

load(meanPFile)
tmpX = v.x;
tmpZ = v.prevFv.zt(end,:);
ii = find(tmpZ==max(tmpZ));
tmpX = tmpX-tmpX(ii);
modP = interp1(newPArcDists,meanMPR,tmpX);
modP = modP./max(modP);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(v.x./1e3,vUplift.u*1e3,'-r','linewidth',2)
xlabel('X (km)')
ylabel('Uplift (mm/yr)')
set(gca,'fontsize',12)
set(gca,'xcolor','k')
set(gca,'ycolor','r')
yy = ylim;
ylim([0,yy(2)]);

subplot(2,2,2)
plot(v.x./1e3,modP,'-b','linewidth',2)
xlabel('X (km)')
ylabel('Scaled Precipitation Rate')
set(gca,'fontsize',12)
set(gca,'xcolor','k')
set(gca,'ycolor','b')
yy = ylim;
ylim([0,yy(2)]);

orient landscape
set(gcf,'PaperType','tabloid')
set(gcf,'Renderer','painters')

saveas(gcf,[saveFol,'Forcing_Functions.fig'],'fig');
saveas(gcf,[saveFol,'Forcing_Functions.pdf'],'pdf');
close

%% Run Analysis
inputs.topoMaxHeight = round(max(useZs(:)),-2)+200;

for i = size(useZs,1):-1:1
    inputs.Topo = useZs(i,:);
    inputs.plotTitle = stepTitles{i};

    sigmas = Calculate_Crustal_Stress_2D(inputs);
    
    save([saveFol,stepTitles{i},'.mat'])

    orient landscape
    set(gcf,'PaperType','tabloid')
    set(gcf,'Renderer','opengl')
    saveas(gcf,[saveFol,stepTitles{i},'_trajectories_bitmap.fig'],'fig');
    print([saveFol,stepTitles{i},'_trajectories_bitmap'],'-dpdf','-r300');
    close

    orient landscape
    set(gcf,'PaperType','tabloid')
    set(gcf,'Renderer','painters')
    saveas(gcf,[saveFol,stepTitles{i},'_trajectories_vector.fig'],'fig');
    saveas(gcf,[saveFol,stepTitles{i},'_trajectories_vector.pdf'],'pdf');
    close

    orient landscape
    set(gcf,'PaperType','tabloid')
    set(gcf,'Renderer','painters')
    saveas(gcf,[saveFol,stepTitles{i},'_trajectories.fig'],'fig');
    saveas(gcf,[saveFol,stepTitles{i},'_trajectories.png'],'png');
    close

    close
    close
end