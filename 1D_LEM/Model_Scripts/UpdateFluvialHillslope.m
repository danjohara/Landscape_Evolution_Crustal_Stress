function [z,eT,eR,basC,L] = UpdateFluvialHillslope(v,er,L,z,basC,boundTypes)
% Name: UpdateFluvialHillslope
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Calculate new model elevation using implicit discretization
% of both the stream power law (fluvial) and diffusion (hillslope) regimes.
% Also update basin statistics.

% Input:
%     v:            Model input structure.
%     er:           Model erosion structure.
%     L:            Variable-lithology structure.
%     z:            Current model elevation.
%     basC:         Current model basin structure.
%     boundTypes:   Grid of model fluvial boundary types.

% Output:
%     z:    New model elevation.
%     eT:   Total eroded material in timestep (excludes landslides).
%     eR:   Timestep erosion rate (excludes landslides).
%     basC: New model basin structure.
%     L:    Variable-lithology structure with updated layer contact
%               elevations.

% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.

% Create implicit scheme coefficients;
% Implicit scheme solves equation b = Az, where b = z0 + u; A = coefficients;
% z = new elevations.
[kappa,psi] = determineKappaPsi(v,er,L,z); % Determine k & D values to use.
alpha = kappa*er.ka^er.m*v.dt/v.dx;
beta = psi*v.dt/v.dx^2;
b = z' + v.u'*v.dt;

% Create Precipitation scaling array.
modX = v.xRange(1):v.dx:v.xRange(2);
if length(er.Ps)==1
    modPrec = ones(size(modX))*er.Ps;
else
    modPrec = interp1(er.Px,er.Ps,modX);
end

% Create array and fill with basin lengths from respective ridge and flags 
% for what each node is: 0 - boundary, 1 - channel, 2 - hillslope, 3 -
% ridge.
flags2 = zeros(size(z)); % 
refX2 = zeros(size(z));
hillslopeBound = zeros(size(z)); %-1: ridge to left, 1: ridge to right.
for i = 1:length(basC.singleBRidgeX)
    if ~isempty(basC.singleBHSNodesI{i}) && isempty(basC.singleBHSNodesI{i})
        if basC.singleBRidgeX(i) < basC.singleBBLX(i)
            hillslopeBound(basC.singleBBLI(i)) = -1;
        else
            hillslopeBound(basC.singleBBLI(i)) = 1;
        end
    end
    
    % Set hillslope locations.
    if ~isempty(basC.singleBHSNodesI{i})
        refX2(basC.singleBHSNodesI{i}) = basC.singleBRidgeX(i) - v.x(basC.singleBHSNodesI{i});
        flags2(basC.singleBHSNodesI{i}) = 2;
    end
    
    % Set ridge locations.
    if ~isempty(basC.singleBRNodesI{i})
        flags2(basC.singleBRNodesI{i}) = 2;
    end
    
    % Set channel locations.
    if ~isempty(basC.singleBCNodesI{i})
        refX2(basC.singleBCNodesI{i}) = basC.singleBRidgeX(i) - v.x(basC.singleBCNodesI{i});
        flags2(basC.singleBCNodesI{i}) = sign(refX2(basC.singleBCNodesI{i}));
        
        % Boundary occurs at channel bottom.
        if sign(refX2(basC.singleBCNodesI{i}(1))) == 1
            flags2(basC.singleBCNodesI{i}(1)) = 0;
        else
            flags2(basC.singleBCNodesI{i}(end)) = 0;
        end
    end
end

% Create coefficient vectors for matrix.
mainDiag = ones(1,length(z));
mainDiag(flags2 == 2) = 1 + 2*beta(flags2 == 2);
% mainDiag(flags2 == 1) = 1 + alpha(flags2==1).*abs(refX2(flags2==1)).^(er.h*er.m);
% mainDiag(flags2 == -1) = 1 + alpha(flags2==-1).*abs(refX2(flags2==-1)).^(er.h*er.m);
mainDiag(flags2 == 1) = 1 + alpha(flags2==1).*abs(refX2(flags2==1).*modPrec(flags2==1)).^(er.h*er.m);
mainDiag(flags2 == -1) = 1 + alpha(flags2==-1).*abs(refX2(flags2==-1).*modPrec(flags2==-1)).^(er.h*er.m);
mainDiag(abs(hillslopeBound)==1) = 1+beta(abs(hillslopeBound)==1);
upDiag = zeros(1,length(z));
downDiag = upDiag;
upDiag(flags2 == 2) = -beta(flags2 == 2);
downDiag(flags2 == 2) = -beta(flags2 == 2);
% upDiag(flags2 == -1) = -alpha(flags2==-1).*abs(refX2(flags2==-1)).^(er.h*er.m);
% downDiag(flags2 == 1) = -alpha(flags2==1).*abs(refX2(flags2==1)).^(er.h*er.m);
upDiag(flags2 == -1) = -alpha(flags2==-1).*abs(refX2(flags2==-1).*modPrec(flags2==-1)).^(er.h*er.m);
downDiag(flags2 == 1) = -alpha(flags2==1).*abs(refX2(flags2==1).*modPrec(flags2==1)).^(er.h*er.m);
upDiag(hillslopeBound==1) = -beta(hillslopeBound==1);
downDiag(hillslopeBound==-1) = -beta(hillslopeBound==-1);

% Where channel boundary corresponds to a boundary flag of 0 (Dirichlet),
% the new elevation is the same as the old elevation
[~,c] = find(flags2==0);
for i = 1:length(c)
    if boundTypes(c(i)) == 0
        b(c(i)) = z(c(i));
    end
end

% Combine coefficient vectors into coefficient matrix for calculation.
A = diag(mainDiag)+diag(upDiag(1:end-1),1)+diag(downDiag(2:end),-1);

% Perform calculation and reorient elevation array.
z0 = z;
z = [A\b]';
z = round(z,6); % This rounding attempts to remove possible numerical errors.

% Calculate erosion total and rate.
eT = z0 + v.u*v.dt - z;
eR = eT/v.dt;

% Update Layer Boundaries
L.zLp = L.zLp + v.u*v.dt;
L.zLn = L.zLn + v.u*v.dt;

% Update Basin Information
basC = UpdateBasins(v,z,basC,eT,eR);
if v.landslide==1
    z = LandslideErode(v,z,basC);
end
end