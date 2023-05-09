function z = LandslideErode(v,z,basC)
% Name: LandslideErode
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Adjust elevations to critical slope, simulating landsliding.

% Input:
%     v:    Model input structure.
%     z:    Current model elevation.
%     basC: Model basin structure.

% Output:
%     z:    New model elevation.

% Start from basin mouths and adjust each upstream node's elevation that
% has a slope greater than critA.
for i = 1:length(basC.singleBRidgeX)
    tmp = abs(v.x-basC.singleBRidgeX(i));
    [~,ridgeI] = find(tmp == min(tmp),1);

    if basC.singleBRidgeX(i) > basC.singleBBLX(i)
        iter = 1;
    elseif basC.singleBRidgeX(i) < basC.singleBBLX(i)
        iter = -1;
    end

    for j = basC.singleBBLI(i)+iter:iter:ridgeI
        alpha = min([v.critA,atand((z(j)-z(j-iter))/v.dx)]);
        z(j) = z(j-iter) + tand(alpha)*v.dx;
    end
end            
end