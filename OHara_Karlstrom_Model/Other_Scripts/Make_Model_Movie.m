%% 
% Name: Make_Model_Movie
% Author: Daniel O'Hara
% Date: 05/09/2023
%
% Description: Script to create the movie for the coupled landscape 
% evolution - crustal stress model for a Cascades-like topography and 
% precipitation shown in O'Hara & Karlstrom (in review).
%
% Reference:
% O'Hara, D., and Karlstrom, L., in review, Distributed volcano erosion 
%   patterns imply coupled volcanism and regional climate in the Cascades 
%   arc, U.S.A. Frontiers in Earth Science.

%% Setup
modFile = '.\LEM_Model_Results\L2_Precip.mat';
filePath = '.\Model_Movie\';
meanPFile = '.\Cascades_Data\Cascades_meanZP.mat';

movName = 'Cascades_Topographic_Stresses';
FR = 5;
movStep = 5;
projStep = 20;

maxTime = 100e6;

phases = {'Initial Topography','Magmatic Construction','Climate Gradient'};
addpath(genpath('.\..\'))

%% Create inputs
load(modFile)
v = v_uplift;
% Elevation, time, phase
allZs = [v.prevFv.zt(:,1:2:end);fv_all.zt(2:end,1:2:end)];
allTs = [v.prevFv.ts;fv_all.ts(2:end)+v.prevFv.ts(end)];
allTs(1) = allTs(2);
allPhases = ones(size(allTs))*2;
allPhases(1) = 1;
allPhases(length(v.prevFv.ts)+1:end) = 3;

% Uplift
tmpUplift1 = v.u(1:2:end);
tmpUplift2 = v.prevFv.v{1}.u(1:2:end);
allUplift = [tmpUplift1;...
    repmat(tmpUplift2,[size(v.prevFv.zt,1)-1,1]);...
    repmat(tmpUplift1,[size(fv_all.zt,1)-1,1])];

% Precipitation
load(meanPFile)
tmpX = v.x;
tmpZ = v.prevFv.zt(end,:);
ii = find(tmpZ==max(tmpZ));
tmpX = tmpX-tmpX(ii);
modP = interp1(newPArcDists,meanMPR,tmpX);
modP = modP./max(modP);

tmpPrecip1 = ones(size(allZs(1,:)));
tmpPrecip2 = modP(1:2:end);
allPrecip = [tmpPrecip1;...
    repmat(tmpPrecip1,[size(v.prevFv.zt,1)-1,1]);...
    repmat(tmpPrecip2,[size(fv_all.zt,1)-1,1])];

x = v.x(1:2:end);
dx = abs(x(2)-x(1));
dz = dx*5;

maxTopos = max(allZs,[],2);
meanTopos = mean(allZs,2);
maxUplift = max(allUplift,[],2);
minPrecip = min(allPrecip,[],2);

%% Create video writer and set up figure window.
vi = VideoWriter([filePath,'\',movName,'.mp4']);
vi.FrameRate = FR;
open(vi);
ff = figure;
set(ff,'Position',get(0,'Screensize'))
colormap(hsv)

%% Create Movie
sPs = {};
stepCounter = 1;
sinceStepCounter = 0;
sPTransparancy = [1,.6,.3,.05];

useX = (x-mean(x))./1e3;

maxDepth = 10100;
minDepth = 100;
tmp = [30e3:5e3:50e3]';
sL = [tmp,ones(size(tmp))*maxDepth-50];

timeStepForcingPos = [0.130732064421669	0.582237355157677	0.334860907759883	0.183661132878385];
timeStepElevationPos = [0.130000000000000	0.409632352941176	0.334659090909091	0.160151819720694];
timeStepStressPos = [0.131771595900439	0.254915095318490	0.333821376281113	0.183661132878385];
cbPos = [0.470473401462991,0.286330935251799,0.015617374525296,0.122520863310793];
allElevationPos = [0.582786004259285	0.576481959474224	0.334659090909092	0.191683508151676];
allForcingPos = [0.582786004259285	0.291238095741680	0.334659090909091	0.183661132878386];

useTime = min([maxTime,allTs(end)]);

for i = 1:movStep:length(allTs)
    disp(sprintf('%d / %d',i,length(allTs)))
    if allTs(i) > maxTime
        break;
    end
    [sigmas,inputs] = Determine_Stress_FromTopo(allZs(i,:),dz,x,dx,minDepth,maxDepth,sL);
    if stepCounter == 1 || sinceStepCounter == projStep
        sinceStepCounter = 0;
        sPs = [{sigmas.sourceProjections};sPs];
        if length(sPs) > length(sPTransparancy)
            sPs(length(sPTransparancy):end) = [];
        end
    end
    
%% Elevations through time
    subplot(3,2,2);
    cla
    hold on
    plot(allTs,maxTopos,'-k','linewidth',2)
    plot(allTs,meanTopos,'--k','linewidth',2)
    xlim([min(allTs),useTime])
    yl = ylim;
    fill([allTs(1),allTs(i),allTs(i),allTs(1),allTs(1)],...
        [yl(1),yl(1),yl(2),yl(2),yl(1)],'k','facealpha',.1,'edgecolor','k','linewidth',.5)
    ylim(yl)
    set(gca,'xscale','log')
    setAxis(gca)
%     xlabel('Time (yr)')
    set(gca,'xtick',[])
    ylabel('Elevation (m)')
    legend('Max','Mean','location','best')
    title('Model Elevations')
    sTopoTime = gca;

%% Precip & Uplift through time
    subplot(3,2,4)
    cla
    yyaxis left
    hold on
    plot(allTs,maxUplift.*1e3,'-r','linewidth',2)
    set(gca,'ycolor','r')
    set(gca,'fontsize',12)
    ylabel('Maximum Uplift (mm/yr)')
    xlim([min(allTs),useTime])
    set(gca,'xscale','log')
    xlim([min(allTs),useTime])

    yyaxis right
    hold on
    plot(allTs,minPrecip,'-b','linewidth',2)
    set(gca,'ycolor','b')
    set(gca,'fontsize',12)
    set(gca,'xcolor','k')
    ylabel('Minimum Precipitation Scaling')
    xlim([min(allTs),max(allTs)])
    xlabel('Time (yr)')
    xlim([min(allTs),useTime])
    
    ylim([0,1])
    yl = ylim;
    fill([allTs(1),allTs(i),allTs(i),allTs(1),allTs(1)],...
        [yl(1),yl(1),yl(2),yl(2),yl(1)],'k','facealpha',.1,'edgecolor','k','linewidth',.5)
    ylim(yl)
    set(gca,'xscale','log')
    box on
    sForcingTime = gca;

    title('Forcing Functions')

%% Forcing Functions at timestep
    subplot(3,2,1)
    cla
    yyaxis left
    plot(useX,allUplift(i,:)*1e3,'-r','linewidth',2);
    set(gca,'ycolor','r')
    set(gca,'fontsize',12)
    ylabel('Uplift (mm/yr)')

    yyaxis right
    plot(useX,allPrecip(i,:),'-b','linewidth',2);
    set(gca,'ycolor','b')
    set(gca,'fontsize',12)
    ylabel('Precipitation Scaling')
    ylim([0,1.1])
    set(gca,'xtick',[]);
    box on
    useT = allTs(i);
    if i == 1
        useT = 1;
    end
    title(sprintf('%s (T = %.3f Myr)',phases{allPhases(i)},useT./1e6))
    sForcingStep = gca;
%% Topography at timestep
    subplot(3,2,3)
    cla
    plot(useX,allZs(i,:)./1e3,'-k','linewidth',2)
    setAxis(gca)
    ylim([0,max(allZs(:))+100]./1e3);
    ylabel('Elevation (km)')
    set(gca,'xtick',[]);
    set(gca,'ytick',[1,2])

    sTopoStep = gca;
%% Crustal Stresses at timestep
    subplot(3,2,5)
    cla
    hold on
    pcolor(sigmas.X./1e3,-sigmas.Z./1e3,sigmas.sigma_1./1e6);shading flat; %axis image
    xlabel('X (km)')
    ylabel('Depth (km)')
    setAxis(gca)
    axis image

    quiver(sigmas.Xcut./1e3,-sigmas.Zcut./1e3,...
        sigmas.sigmaP_dx_cut./1e3,sigmas.sigmaP_dz_cut./1e3,0,...
        'color',[1,1,1]*.5,'marker','none','ShowArrowHead','off','linewidth',inputs.VectorWidth)

    for j = 1:length(sPs)
        for k = 1:length(sPs{j})
            p1 = patchline(sPs{j}{k}(:,1)./1e3,-sPs{j}{k}(:,2)./1e3,...
                'edgecolor','k','linewidth',2,'edgealpha',sPTransparancy(j));
        end
    end

    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')

    set(gca,'position',timeStepStressPos)
    set(cb,'position',cbPos)
    
    caxis(inputs.magScale2)

%% Set Positions
    set(sTopoTime,'position',allElevationPos);
    set(sForcingTime,'position',allForcingPos);
    set(sForcingStep,'position',timeStepForcingPos)
    set(sTopoStep,'position',timeStepElevationPos)

    %% Wrap up
    stepCounter = stepCounter+1;
    sinceStepCounter = sinceStepCounter+1;

    frame = getframe(ff);
    writeVideo(vi,frame)
end

% Close writer.
close(vi)
close

%% 
function setAxis(gg)
    set(gg,'xcolor','k')
    set(gg,'ycolor','k')
    set(gg,'fontsize',12)
    set(gg,'box','on')
end