function [sigmaXX,sigmaYY,sigmaZZ,sigmaXY,sigmaXZ,sigmaYZ] = Boussonesq_Convolution(x,y,z,Fv,topoLoad,lambda,G,kernel_radius,kernel_res)
% Name: Boussonesq_Convolution
% Original Python Authors: Richard H. Styron and Eric A. Hetland
% Matlab Adaptation: Daniel O'Hara
% Date: 11/20/2022
%
% Description: Matlab recreation of the Halfspace Python module from Styron
% & Hetland (2015)(https://github.com/cossatot/halfspace). Script performs 
% the Boussonesq convolution for topographic loading. Used in a coupled 
% landscape evolution - crustal stress model by O'Hara & Karlstrom (in 
% review).
%
% Input:
%   x:              Grid x-coordinates.
%   y:              Grid y-coordinates.
%   z:              Grid z-coordinates.
%   Fv:             Gravitational force.
%   topoLoad:       Topography.
%   lambda:         First Lame parameter.
%   G:              Second Lame parameter.
%   kernel_radius:  Kernel radius for convolution.
%   kernel_res:     Kernel resolution (grid resolution).
%
% Ouput:
%   sigmaXX:        Sigma_x in x-direction.
%   sigmaYY:        Sigma_y in y-direction.
%   sigmaZZ:        Sigma_z in z-direction.
%   sigmaXY:        Sigma_x in y-direction.
%   sigmaXZ:        Sigma_x in z-direction.
%   sigmaYZ:        Sigma_y in z-direction.
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

stressList = {'XX','YY','ZZ','XY','XZ','YZ'};
if isnan(y)
    sigmaXX = zeros(length(z),length(x));
    sigmaZZ = sigmaXX;
    sigmaXZ = sigmaXX;

    sigmaYY = NaN;
    sigmaXY = NaN;
    sigmaYZ = NaN;
else
    sigmaXX = zeros(length(y),length(x),length(z));
    sigmaYY = sigmaXX;
    sigmaZZ = sigmaXX;
    sigmaXY = sigmaXX;
    sigmaXZ = sigmaXX;
    sigmaYZ = sigmaXX;
end

% Run through each component and depth.
for i = 1:length(stressList)
    if isnan(y) && (i==2 || i==4 || i==6)
        continue;
    end

    for j = 1:length(z)
        %% Make Kernel
        kernel = Make_Boussonesq_Kernel(y,z(j),stressList{i},Fv,lambda,G,kernel_radius,kernel_res);
        
        %% Perform Convolution
        if isnan(y)
            tmp = conv(topoLoad,kernel,'same');
            evalc(sprintf('sigma%s(j,:) = tmp;',stressList{i}));
        else
            tmp = conv2(topoLoad,kernel,[],'same');
            evalc(sprintf('sigma%s(:,:,j) = tmp;',stressList{i}));
        end
    end
end