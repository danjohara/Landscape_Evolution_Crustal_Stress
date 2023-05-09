function [sigma1,sigma3,sigmaP] = Calculate_Principal_Stresses_2D(sigmaXX,sigmaZZ,sigmaXZ)
% Name: Calculate_Principal_Stresses_2D
% Author: Daniel O'Hara
% Date: 11/20/2022
%
% Description: Script to calculate principal stresses for an x-z profile.
%
% Input:
%   sigmaXX:        Sigma_x in x-direction.
%   sigmaZZ:        Sigma_z in z-direction.
%   sigmaXZ:        Sigma_x in z-direction.
%
% Output:
%   sigma1:        Sigma_1 principal stress.
%   sigma3:        Sigma_3 principal stress.
%   sigmaP:        Sigma_1 stress orientation.

sigmaP = atand(2*sigmaXZ./(sigmaXX-sigmaZZ))./2;
sigma1 = (sigmaXX+sigmaZZ)./2 + sqrt(((sigmaXX-sigmaZZ)./2).^2 + sigmaXZ.^2);
sigma3 = (sigmaXX+sigmaZZ)./2 - sqrt(((sigmaXX-sigmaZZ)./2).^2 + sigmaXZ.^2);
end