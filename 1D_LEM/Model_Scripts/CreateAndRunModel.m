function fv_all = CreateAndRunModel(f,v,en,er,P,L,U)
% Name: CreateAndRunModel
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Uses the structured model variables (f,v,en,er,P,L) to run
% the 1D landscape evolution model from O'Hara et al. (2019), which allows 
% for different uplift, precipitation lithology parameters. Afterwards,  
% saves the results to the given file.

% Input Model Structures:
%   f: Save file parameters.
%       filepath:   Save file path.
%       filename:   Name of the save file.
%   v: Main model parameters.
%       dx:         Model grid spacing.
%       nx:         Number of spatial grid nodes.
%       xRange:     Spatial node range.
%       allPlot:    Flag to plot model steps:
%           0:          Plot all landscape with no distinction.
%           1:          Plot only channels with distinction.
%           2:          Plot entire landscape with channel distinction.
%       savePlot:   Flag to save plots.
%       sinkHeight: Base level height.
%       landslide:  Flag to enable landslides, ensuring landscape does not 
%                       reach slopes above critA.
%       critA:      Critical slope angle (in degrees) for landsliding.
%       xr0:        Initial ridge location (assuming steady state).
%       minTime:    Minimum model run time.
%       ssRHeightPer: Model stop criteria for how close model divide height 
%                       is to the initial divide.
%       ssRPosPer:  Model stop criteria for how close model divide position
%                       is to the initial divide.
%       initRHeight: Initial divide height for stop criteria.
%       initRPos:   Initial divide position for stop criteria.
%   en: End model parameters.
%       st:         Model output save timestep.
%       tf:         Model final time.
%       dt:         Model timestep.
%   er: Erosion model parameters.
%       m:          Stream power drainage area dependence exponent. 
%       n:          Stream power slope dependence exponent (model setup
%                       requires n = 1).
%       h:          Hack's law exponent (~1.67).
%       ka:         Hack's law coefficient (~.57).
%       u0:         Background uplift rate.
%       k:          Steam power erosion coefficient.
%       D:          Hillslope diffusion coefficient (model setup uses
%                       linear diffusion).
%       Ps:         Precipitation scaling array.
%       Px:         Precipitation scaling coordinates.
%   P: Uplift perturbation model parameters, using uplift fuction from 
%           O'Hara et al. (2019).
%       runPerturbationModel:  Flag for whether uplift perturbation is
%                       modeled.
%       xP:         Location of perturbation maximum uplift.
%       R:          Perturbation radius.
%       umax:       Perturbation maximum uplift rate.
%       b:          Perturbation uplift function exponent.
%       uT:         Perturbation uplift timescale.
%       dt:         Perturbation uplift timestep.
%       st:         Model output save timestep during perturbation uplift.
%       maxuT:      Maximum perturbation uplift (for syncing of multiple
%                       perturbation models).
%   L: Varying lithology model parameters.
%       runLithologyModel:  Flag for whether lithology changes is
%                       modeled.
%       plotLithology:  Flag for whether to plot different layer.
%       kL:             Layer erodibility value (for channels).
%       DL:             Layer diffusvity value (for hillslopes).
%       wL:          	Layer width.
%       xL:             Spatial coordinate of the layer's center position.
%       thetaL:         Layer dipping angle (model does not allow 90
%                           degrees).
%       dL:             Initial depth of the layer's center position
%                           (assuming postive downward).
%       Thickness:      Layer half thickness (assumes constant value).
%   U: Uplift gradient model parameters.
%       runUpliftGradientModel: Flag for whether uplift gradient is
%                               modeled.
%       leftU:  Uplift rate of left model boundary (assuming linear
%           gradient).
%       rightU: Uplift rate of right model boundary (assuming linear
%           gradient).
% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.
    
%% Setup Model
f.filepath = strrep(f.filepath,'\','/'); % Make sure filepath fits script notation.
% Create Initial Landscape:
[v.x,v.l,v.z,v.basC,v.lc] = CreateInitialProfile(v,er);
% v.x: model grid.
% v.l: model grid distances from ridge.
% v.z: model elevations.
% v.basC: model basin variables.
% v.lc: critical length for fluvial regime.

%% Create Uplift Field
v = createUpliftField(v,P,er,U);
v0 = v;

%% Initialize Lithology Boundaries
L = CreateInitialLithologyLayers(v,L);

%% Run model
% If perturbation model is being run, perform three consecutive models that
% 1) Uplift the perturbation, 2) Uplift landscape with the background
% uplift using the same timestep as the perturbation (for syncing across
% multiple models), and 3) Uplift landscape with background uplift until
% final time is reached. Otherwise, run just one model.
if P.runPerturbationModel
    %% Perturbation Uplift
    v.tf = P.uT;
    v.st = P.st;
    v.dt = P.dt;
    v_uplift = v;
    L_uplift = L;
    fv_uplift = RunModel(v,er,L);
    
    v.z = fv_uplift.zt(end,:);
    v.basC = fv_uplift.basCs{end};
    v.u = v.u0;
    
    L.zLp = fv_uplift.zLps{end};
    L.zLn = fv_uplift.zLns{end};
    
    %% Finish Perturbation Uplift Time
    if P.uT ~= P.maxuT
        v.tf = P.maxuT - P.uT;
        v.minTime = v.minTime - P.uT;
        v_tmp = v;
        L_tmp = L;
        fv_tmp = RunModel(v,er,L);
        
        v.z = fv_tmp.zt(end,:);
        v.basC = fv_tmp.basCs{end};
        
        L.zLp = fv_tmp.zLps{end};
        L.zLn = fv_tmp.zLns{end};
    else
        fv_tmp = [];
        v_tmp = [];
        L_tmp = [];
    end
    
    %% Run Model to End Time
    v.tf = en.tf;
    v.dt = en.dt;
    v.st = en.st;
    v.minTime = v0.minTime - P.maxuT;
    v_erode = v;
    L_erode = L;
    fv_erode = RunModel(v,er,L);
    
    %% Combine FV's
    if isempty(fv_tmp)
        allZT = [fv_uplift.zt;fv_erode.zt(2:end,:)];
        allLT = [fv_uplift.lt;fv_erode.lt(2:end,:)];
        allET = [fv_uplift.et;fv_erode.et(2:end,:)];
        allER = [fv_uplift.er;fv_erode.er(2:end,:)];
        allTs = [fv_uplift.ts;fv_erode.ts(2:end)+fv_uplift.ts(end)];
        allBCs = [fv_uplift.basCs;fv_erode.basCs(2:end)];
        allZLps = [fv_uplift.zLps;fv_erode.zLps(2:end)];
        allZLns = [fv_uplift.zLns;fv_erode.zLns(2:end)];
        allVs = {v_uplift;v_erode};
    else
        allZT = [fv_uplift.zt;fv_tmp.zt(2:end,:);fv_erode.zt(2:end,:)];
        allLT = [fv_uplift.lt;fv_tmp.lt(2:end,:);fv_erode.lt(2:end,:)];
        allET = [fv_uplift.et;fv_tmp.et(2:end,:);fv_erode.et(2:end,:)];
        allER = [fv_uplift.er;fv_tmp.er(2:end,:);fv_erode.er(2:end,:)];
        allTs = [fv_uplift.ts;fv_tmp.ts(2:end)+fv_uplift.ts(end);fv_erode.ts(2:end)+fv_uplift.ts(end)+fv_tmp.ts(end)];
        allBCs = [fv_uplift.basCs;fv_tmp.basCs(2:end);fv_erode.basCs(2:end)];
        allZLps = [fv_uplift.zLps;fv_tmp.zLps(2:end);fv_erode.zLps(2:end)];
        allZLns = [fv_uplift.zLns;fv_tmp.zLns(2:end);fv_erode.zLns(2:end)];
        allVs = {v_uplift;v_tmp;v_erode};
    end
    
    fv_all.zt = allZT;
    fv_all.lt = allLT;
    fv_all.et = allET;
    fv_all.er = allER;
    fv_all.ts = allTs;
    fv_all.basCs = allBCs;
    fv_all.zLps = allZLps;
    fv_all.zLns = allZLns;
else
    %% Run single model
    v.st = en.st;
    v.tf = en.tf;
    v.st = en.st;
    v.dt = en.dt;
    
    fv_all = RunModel(v,er,L);
    allVs = {v};
end
fv_all.v = allVs;

%% Combine all model output structures into one
if P.runPerturbationModel
    save([f.filepath,'/',f.filename,'.mat'],'v_erode','v_tmp','v_uplift','fv_all','P','L','en','er','f');
else
    save([f.filepath,'/',f.filename,'.mat'],'v','fv_all','P','L','en','er','f');
end
end