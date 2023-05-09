function [kappa,psi] = determineKappaPsi(v,er,L,z)
% Name: determineKappaPsi
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Determine the channel erobility (k) and hillslope 
%   diffusivity (D) landscape values, for the case of variable lithology.

% Input:
%     v:    Model input structure.
%     er:   Model erosion structure.
%     L:    Variable-lithology structure.
%     z:    Current elevation grid.

% Output:
%     kappa: Channel erodibility grid.
%     psi:   Hillslope diffusivity grid.

    % Set grid values to background k & D.
    kappa = ones(size(v.x))*er.k;
    psi = ones(size(v.x))*er.D;
    if ~L.runLithologyModel
        return
    end
    
    % Find locations where the landscape surface intersects with the layer.
    ti1 = find(~isnan(L.zLp),1);
    ti2 = find(~isnan(L.zLp),1,'last');
    ti3 = find(~isnan(L.zLn),1,'last');
    ti4 = find(~isnan(L.zLn),1);
    lX = [v.x(ti1:ti2),fliplr(v.x(ti4:ti3)),v.x(ti1)];
    lZ = [L.zLp(ti1:ti2),fliplr(L.zLn(ti4:ti3)),L.zLp(ti1)];
    
    inP = inpolygon(v.x,z,lX,lZ);
    
%     t1 = z<=L.zLp;
%     t2 = z>=L.zLn;
%     t3 = t1.*t2;
    kappa(inP==1) = L.kL;
    psi(inP==1) = L.DL;
    
end