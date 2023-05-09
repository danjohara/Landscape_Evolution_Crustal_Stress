function L = CreateInitialLithologyLayers(v,L)
% Name: CreateInitialLithologyLayers
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Determines the initial upper and lower contacts of a 
%   rock layer of different erodibility within an otherwise-uniform 
%   lithology. Vertical layers are currently not allowed (theta != 90).

% Input:
%     v:    Model input structure.
%     L:    Variable-lithology structure.

% Output:
%     fv:   Structured output variable.

L.zLp = NaN;
L.zLn = NaN;
if ~L.runLithologyModel
    return
else
    % Determine initial landscape elevation at xL (symmetry point of
    % layer)
    z0_xL = interp1(v.x,v.z,L.xL);
    
    if L.thetaL ~= 90 && L.thetaL ~= -90
        % Determine local horizontal coordinates of global grid for upper 
        % and lower contacts.
        xLp = (v.x - L.xL + sind(L.thetaL)*L.Thickness)./cosd(L.thetaL);
        xLn = (v.x - L.xL + sind(L.thetaL)*-L.Thickness)./cosd(L.thetaL);
        
        % Using local horizontal coordinates, calculate contact elevations
        % in global domain.
        L.zLp = z0_xL - L.dL + sind(L.thetaL)*xLp + cosd(L.thetaL)*L.Thickness;
        L.zLn = z0_xL - L.dL + sind(L.thetaL)*xLn + cosd(L.thetaL)*-L.Thickness;
        
        % Set everything outside the bounds to NaN.
        L.zLp(abs(xLp) >= L.wL/2) = NaN;
        L.zLn(abs(xLn) >= L.wL/2) = NaN;
        
    else
        % Determine elevation of the upper and lower contact points.
        L.zLp = ones(size(v.x))*(z0_xL - L.dL + L.wL/2);
        L.zLn = ones(size(v.x))*(z0_xL - L.dL - L.wL/2);
        
        % Set everything outside the bounds to NaN.
        L.zLp(abs(v.x-L.xL)>L.Thickness) = NaN;
        L.zLn(abs(v.x-L.xL)>L.Thickness) = NaN;
    end
end
end