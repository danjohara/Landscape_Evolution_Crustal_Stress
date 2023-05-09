function [sigmaP_dx,sigmaP_dz] = Get_Stress_Orientations(sigmaXX,sigmaZZ,sigmaXZ,sigmaP_scale)
% Name: Get_Stress_Orientations
% Author: Daniel O'Hara
% Date: 11/20/2022
%
% Description: Calculate sigma_1 stress orientations as a vector of x-z
%   data. Used for Matlab's quiver function.
%
% Input:
%   sigmaXX:        Sigma_x in x-direction.
%   sigmaZZ:        Sigma_z in z-direction.
%   sigmaXZ:        Sigma_x in z-direction.
%   sigmaP_scale:   Vector scaling value for stress orientation plots.
%
% Input:
%   sigmaP_dx:        Sigma_1 principal stress vector x-value.
%   sigmaP_dz:        Sigma_3 principal stress vector x-value.

[~,~,sigmaP] = Calculate_Principal_Stresses_2D(sigmaXX,sigmaZZ,sigmaXZ);

sigmaP_dx = zeros(size(sigmaP))*NaN;
sigmaP_dz = sigmaP_dx;

t1 = sigmaP < 0;
t2 = sigmaP > 0;
t3 = sigmaXX-sigmaZZ<0;
t4 = sigmaXX-sigmaZZ>0;

tt1 = t1.*t4;
tt2 = t2.*t4;

sigmaP_dx(tt1==1) = sind(90-abs(sigmaP(tt1==1)))*sigmaP_scale;
sigmaP_dx(tt2==1) = -sind(90-abs(sigmaP(tt2==1)))*sigmaP_scale;
sigmaP_dz(t4) = cosd(90-abs(sigmaP(t4)))*sigmaP_scale;

tt1 = t1.*t3;
tt2 = t2.*t3;

sigmaP_dx(tt1==1) = -sind(abs(sigmaP(tt1==1)))*sigmaP_scale;
sigmaP_dx(tt2==1) = sind(abs(sigmaP(tt2==1)))*sigmaP_scale;
sigmaP_dz(t3) = cosd(abs(sigmaP(t3)))*sigmaP_scale;
end