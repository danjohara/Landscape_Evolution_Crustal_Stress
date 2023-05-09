function sigmas = Calculate_Crustal_Stress_2D(inputs)
% Name: Calculate_Crustal_Stress_2D
% Python Code Authors: Richard H. Styron and Eric A. Hetland
% Matlab Adaptation: Daniel O'Hara
% Date: 11/20/2022
%
% Description: Matlab recreation of the Halfspace Python module from Styron
% & Hetland (2015)(https://github.com/cossatot/halfspace). Script  
% calculates the first-order crustal stresses for topographic loading. Used 
% in a coupled landscape evolution - crustal stress model by O'Hara & 
% Karlstrom (in review).

% Input Model Structures:
%   inputs: Model parameters.
%       dx:             Model horizontal grid resolution.
%       x:              Model horizontal coordinates.
%       dz:             Model vertical grid resolution.
%       maxDepth:       Maximum model vertical grid depth.
%       minDepth:       Minimum model vertical grid depth.
%       Topo:           Model topography grid.
%
%       Analysis:       Analysis flag. 1 = full grid; 2 = point load from 
%                           peak; 3 = point load from centroid (eventually)
%       lambda:         First Lame parameter.
%       G:              Second Lame parameter.
%       rho:            Crust density
%       kernel_radius:  Kernel radius for convolution.
%
%       makePlots:      Flag to make plots.
%       sigmaP_scale:   Vector scaling value for stress orientation plots.
%       VectorWidth:    Vector widths for stress orientation plots.
%       SigmaP_cutByX:  Vector spacing in x-direction for stress
%                           orientation plots.
%       SigmaP_cutByZ:  Vector spacing in z-direction for stress
%                           orientation plots.
%       magScale1:      Color axis for \sigma_x and \sigma_z plots.
%       magScale2:      Color axis for \sigma_1 plots.
%       sourceLocs:     Source x-locations for trajectory plots.
%       topoMaxHeight:  Maximum topographic height for plotting.
%       plotTitle:      Plot titles.
%
% Output:
%   sigmas:
%       X:                  X-grid coordinates.
%       Z:                  Z-grid coordinates.
%       sigma_XX:           Grid of sigma_XX values.
%       sigma_ZZ:           Grid of sigma_ZZ values.
%       sigma_XZ:           Grid of sigma_XZ values.
%       sigma_1:            Grid of sigma_1 values.
%       sigma_3:            Grid of sigma_3 values.
%       sigma_P:            Grid of sigma_1 stress orientations.
%       Xcut:               Cut grid of x-coordinates.
%       Zcut:               Cut grid of y-coordinates.
%       sigmaP_dx_cut:      Grid of sigma_1 x-direction vector componenets.
%       sigmaP_dz_cut:      Grid of sigma_1 z-direction vector componenets.
%       sourceProjections:  x-z sigma-1 trajectory coordinates.
%
% References:
% Styron, R.H., and Hetland, E.A., 2015, The weight of the mountains: 
%   Constraints on tectonic stress, friction, and fluid pressure in the 
%   2008 Wenchuan earthquake from estimates of topographic loading: Journal 
%   of Geophysical Research: Solid Earth, v. 120, p. 2697–2716, 
%   doi:10.1002/2014JB011338
%
% O'Hara, D., and Karlstrom, L., in review, Distributed volcano erosion 
%   patterns imply coupled volcanism and regional climate in the Cascades 
%   arc, U.S.A. Frontiers in Earth Science.
    
%% Setup
g = 9.8;

% Get User Material Properties
G = inputs.G;
lambda = inputs.lambda;
rho = inputs.rho;

% Get User array data
topoLoad = inputs.Topo;
dz = inputs.dz;
dx = inputs.dx;
z = inputs.minDepth:dz:inputs.maxDepth;
x = inputs.x;

% Get User Kernel Data
kernel_rad = inputs.kernel_radius;

% Get User Plotting Values
makePlots = inputs.makePlots;
sigmaP_scale = inputs.sigmaP_scale;
cutByX = inputs.sigmaP_cutByX;
cutByZ = inputs.sigmaP_cutByZ;
magRange1 = inputs.magScale1;
magRange2 = inputs.magScale2;
vectWidth = inputs.VectorWidth;
plotTitle = inputs.plotTitle;
topoMaxHeight = inputs.topoMaxHeight;

% Get Type of Analysis - 1 = full grid; 2 = point load from peak; 3 = point
% load from centroid (eventually)
analysis_Flag = inputs.Analysis;
if analysis_Flag ==1
    topoLoad = -topoLoad;
    pointPos = NaN;
elseif analysis_Flag == 2
    tmpLoad = zeros(size(topoLoad));
    ii = find(topoLoad == max(topoLoad),1);
    tmpLoad(ii) = max(topoLoad);

    topoLoad = -tmpLoad;
    pointPos = x(ii);
end


% Create Grids and Forces
Fv = -rho*g;
[X,Z] = meshgrid(x,z);

% Create Kernel
kernel_len = round(kernel_rad*2/(dx+1));

%% Perform Boussonesq Convolution
[sigmaXX_B,sigmaYY_B,sigmaZZ_B,sigmaXY_B,sigmaXZ_B,sigmaYZ_B] = Boussonesq_Convolution(x,NaN,z,Fv,topoLoad,lambda,G,kernel_rad,dx);

% For now, until the Cerruti correction is included...
sigmaXX = sigmaXX_B;
sigmaYY = sigmaYY_B;
sigmaZZ = sigmaZZ_B;
sigmaXY = sigmaXY_B;
sigmaXZ = sigmaXZ_B;
sigmaYZ = sigmaYZ_B;

%% Perform Cerruti Correctionb (not yet implemented)
% [sigmaXX_C,sigmaYY_C,sigmaZZ_C,sigmaXY_C,sigmaXZ_C,sigmaYZ_C] = Cerruti_Convolution(x,NaN,z,Fv,topoLoad,lambda,G,kernel_rad,dx,sigmaXX_B,sigmaYY_B,sigmaXY_B);

%% Get Principal Stress and Directions
[sigma1,sigma3,sigmaP] = Calculate_Principal_Stresses_2D(sigmaXX,sigmaZZ,sigmaXZ);

%% Determine sigma_1 Vector Orientations
[sigmaP_dx,sigmaP_dz] = Get_Stress_Orientations(sigmaXX,sigmaZZ,sigmaXZ,sigmaP_scale);

Xcut = X(cutByZ:cutByZ:end,cutByX:cutByX:end);
Zcut = Z(cutByZ:cutByZ:end,cutByX:cutByX:end);
sigmaP_dx_cut = sigmaP_dx(cutByZ:cutByZ:end,cutByX:cutByX:end);
sigmaP_dz_cut = sigmaP_dz(cutByZ:cutByZ:end,cutByX:cutByX:end);

%% Project stress directions from source
sourceProjections = {};
if isfield(inputs,'sourceLocs')
    px = inputs.minDepth;

    scatXX = scatteredInterpolant(X(:),Z(:),sigmaXX(:));
    scatZZ = scatteredInterpolant(X(:),Z(:),sigmaZZ(:));
    scatXZ = scatteredInterpolant(X(:),Z(:),sigmaXZ(:));

    for i = 1:size(inputs.sourceLocs,1)
        tmpP = inputs.sourceLocs(i,:);
        while tmpP(end,2) > -px+10
            tmp_sigX = scatXX(tmpP(end,1),tmpP(end,2));
            tmp_sigZ = scatZZ(tmpP(end,1),tmpP(end,2));
            tmp_sigXZ = scatXZ(tmpP(end,1),tmpP(end,2));

            [tmpP_dx,tmpP_dz] = Get_Stress_Orientations(tmp_sigX,tmp_sigZ,tmp_sigXZ,px);
            tmpP = [tmpP;tmpP(end,1)+tmpP_dx,tmpP(end,2)-tmpP_dz];
        end
        sourceProjections = [sourceProjections;{tmpP}];
    end
end

%% Plots
if makePlots
    %% Plot sigma_XX and sigma_ZZ
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    pcolor(X./1e3,-Z./1e3,sigmaXX./1e6); shading flat; %axis imagesigma1
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{XX} (MPA)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    if ~isnan(magRange1)
        caxis(magRange1);
    end
    colormap(bluewhitered)
    axis image
    
    subplot(1,2,2)
    pcolor(X./1e3,-Z./1e3,sigmaZZ./1e6); shading flat; %axis image;
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{ZZ} (MPA)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    if ~isnan(magRange1)
        caxis(magRange1);
    end
    axis image
    
    %% Plot sigma_1 and topography
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    pcolor(X./1e3,-Z./1e3,sigma1./1e6);shading flat; %axis image
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    colormap(hsv)
    axis image
    if ~isnan(magRange2)
        caxis(magRange2);
    end
    
    subplot(1,2,2)
    plot(x./1e3,inputs.Topo./1e3,'-k','linewidth',2)
    xlabel('X (km)')
    ylabel('Topography (km)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    ylim([0,topoMaxHeight./1e3])
    
    %% Plot Orientations
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(3,2,[2,3]*2)
    hold on
    pcolor(X./1e3,-Z./1e3,sigma1./1e6); shading flat; %axis image; 
    quiver(Xcut./1e3,-Zcut./1e3,sigmaP_dx_cut./1e3,sigmaP_dz_cut./1e3,0,'k','marker','none','ShowArrowHead','off','linewidth',vectWidth)
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    if ~isnan(magRange2)
        caxis(magRange2);
    end
    colormap(hsv)
    box on
    axis image
    xx = xlim;
    yy = ylim;
    title(plotTitle)
    
    subplot(3,2,1*2)
    plot(x./1e3,inputs.Topo./1e3,'-k','linewidth',2)
    xlabel('X (km)')
    ylabel('Topography (km)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    box on
    ylim([0,topoMaxHeight./1e3])
    
    subplot(3,2,[2,3]*2-1)
    pcolor(X./1e3,-Z./1e3,sigma1./1e6); shading flat; %axis image; 
    hold on
    for i = 1:length(sourceProjections)
        plot(sourceProjections{i}(:,1)./1e3,-sourceProjections{i}(:,2)./1e3,'-k','linewidth',2)
    end
    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')
    xlabel('X (km)')
    ylabel('Depth (km)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    if ~isnan(magRange2)
        caxis(magRange2);
    end
    colormap(hsv)
    cc = caxis;
    box on
    axis image
    xlim(xx)
    ylim(yy)
    title(plotTitle)
    
    %% Plot Orientations - for vector graphics
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(3,2,[2,3]*2)
    hold on
    caxis(cc)
    quiver(Xcut./1e3,-Zcut./1e3,sigmaP_dx_cut./1e3,sigmaP_dz_cut./1e3,0,'k','marker','none','ShowArrowHead','off','linewidth',vectWidth)
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')
    colormap(hsv)
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    axis image
    box on
    xlim(xx)
    ylim(yy)
    title(plotTitle)
    
    subplot(3,2,1*2)
    plot(x./1e3,inputs.Topo./1e3,'-k','linewidth',2)
    xlabel('X (km)')
    ylabel('Topography (km)')
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    box on
    ylim([0,topoMaxHeight./1e3])
    
    subplot(3,2,[2,3]*2-1)
    hold on
    caxis(cc)
    for i = 1:length(sourceProjections)
        plot(sourceProjections{i}(:,1)./1e3,-sourceProjections{i}(:,2)./1e3,'-k','linewidth',2)
    end
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')
    colormap(hsv)
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    axis image
    box on
    xlim(xx)
    ylim(yy)
    title(plotTitle)
    
    %% Plot background - for bitmap graphic
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(3,2,[2,3]*2)
    hold on
    pcolor(X./1e3,-Z./1e3,(sigma1./1e6)); shading flat; %axis image; 
    caxis(cc)
    xlabel('X (km)')
    ylabel('Depth (km)')
    cb = colorbar;
    ylabel(cb,'\sigma_{1} (MPA)')
    colormap(hsv)
    axis image
    box on
    set(gca,'xcolor','k')
    set(gca,'ycolor','k')
    set(gca,'fontsize',12)
    xlim(xx)
    ylim(yy)
    title(plotTitle)
    if ~isnan(magRange2)
        caxis(magRange2);
    end
end

%% Make Output
sigmas.X = X;
sigmas.Z = Z;
sigmas.sigma_XX = sigmaXX;
sigmas.sigma_ZZ = sigmaZZ;
sigmas.sigma_XZ = sigmaXZ;
sigmas.sigma_1 = sigma1;
sigmas.sigma_3 = sigma3;
sigmas.sigma_P = sigmaP;

sigmas.Xcut = Xcut;
sigmas.Zcut = Zcut;
sigmas.sigmaP_dx_cut = sigmaP_dx_cut;
sigmas.sigmaP_dz_cut = sigmaP_dz_cut;

sigmas.sourceProjections = sourceProjections;
end