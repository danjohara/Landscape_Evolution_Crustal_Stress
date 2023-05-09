function [fv] = RunModel(v,er,L)
% Name: RunModel
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Runs the 1D landscape evolution model using the parameters
% given in the structured input variables (v,er,L).

% Input:
%     v:    Model input structure.
%     er:   Model erosion structure.
%     L:    Variable-lithology structure.

% Output:
%     fv:   Structured output variable.

% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.

% Create uplift grid, if not already done.
if length(v.u) == 1
    v.u = ones(size(v.z))*v.u;
end

% If no timestep given, calculate it from grid spacing using a Courant
% number of 1.
if ~isfield(v,'dt')
    v.dt = v.dx/(v.k*v.ka^v.m*max(v.x));
end

% Ensure timestep is not larger than the saving step and final time.
if v.dt > v.st
    v.dt = v.st;
end

if v.dt > v.tf
    v.dt = v.tf;
end

% Unwrap changing variable structure.
z = v.z;
basC = v.basC;
L0 = L;

% Define boundary types for entire grid. 0 = Dirichlet (base level), 1 =
% Neumann. Dirichlet assumed at ends and Neumann (when applicable)
% everywhere else.
bTypes = ones(size(z));
bTypes(1) = 0;
bTypes(end) = 0;

% Define grid of times that will save model output and display model (when
% plotting).
pt = 0:v.st:v.tf;
ptc = 1;
if pt(end)~= v.tf
    pt = [pt,v.tf];
    spt = length(pt);
else
    spt = length(pt);
end
ptf = zeros(1,spt);

% Initialize save variables.
zt = zeros(spt,length(z));
eRt = zt;
eTt = zt;
lt = zt;
ts = zeros(spt,1);
basCs = cell(spt,1);
zLps = cell(spt,1);
zLns = cell(spt,1);

zt(1,:) = z;
eRt(1,:) = v.u;
lt(1,:) = v.l;
ts(1) = 0;
basCs{1} = basC;
zLps{1} = L.zLp;
zLns{1} = L.zLn;

% Model loop.
i = v.dt;
while i <=v.tf
    disp(i)
    % Update elevations and basin structure.
    [z,eT,eR,basC,L] = UpdateFluvialHillslope(v,er,L,z,basC,bTypes); 

     % Plot model landscape (if applicable).
    if v.allPlot > 0
        PlotXZ(v,L,z,basC,i)
    end
    
    % Check and save data (if applicable).
    if ptc <= length(ptf) && i >= pt(ptc) && ptf(ptc) == 0
        ptf(ptc) = 1;
        ptc = ptc+1;
    
        zt(ptc,:) = z;
        eTt(ptc,:) = eT;
        eRt(ptc,:) = eR;
        lt(ptc,:) = abs(v.x-basC.singleBRidgeX(1));
        ts(ptc) = i;
        basCs{ptc} = basC;
        zLps{ptc} = L.zLp;
        zLns{ptc} = L.zLn;

        % Check and save model plot.
        if v.savePlot == 1
            PlotXZ(v,L,z,basC,i)
            orient landscape
            set(gcf,'PaperType','tabloid')
            saveas(gcf,[v.filepath,'/',v.filename,'.fig']);
        end
    end
    
    % Stop model if landscape state is within of percentage of steady
    % state.
    if i > v.minTime && length(basC.singleBRidgeX) == 2
        if abs(basC.singleBRidgeX(end)-v.initRPos)/v.initRPos <= v.ssRPosPer && abs(basC.singleBRidgeZ(end)-v.initRHeight)/v.initRHeight <= v.ssRHeightPer
            ptc = ptc+1;
            zt(ptc,:) = z;
            eTt(ptc,:) = eT;
            eRt(ptc,:) = eR;
            lt(ptc,:) = abs(v.x-basC.singleBRidgeX(1));
            ts(ptc) = i;
            basCs{ptc} = basC;
            zLps{ptc} = L.zLp;
            zLns{ptc} = L.zLn;
            
            ptc = ptc+1;
            if ptc < length(ts)
                zt(ptc:end,:) = [];
                eRt(ptc:end,:) = [];
                eTt(ptc:end,:) = [];
                lt(ptc:end,:) = [];
                ts(ptc:end) = [];
                basCs(ptc:end) = [];
                zLps(ptc:end) = [];
                zLns(ptc:end) = [];
            end
            break
        end
    end
            

    % Update loop counter, ensure loop stops at Tf.
    if v.tf-i < v.dt && v.tf-i > 0
        i = v.tf;
    else
        i = i+v.dt;
    end    
end

% Create output structure variable.
if isempty(basCs{end})
    basCs = basCs(1:end-1);
end
fv.zt = zt;
fv.et = eTt;
fv.er = eRt;
fv.lt = lt;
fv.ts = ts;
fv.basCs = basCs;
fv.v = v;
fv.zLps = zLps;
fv.zLns = zLns;
end

