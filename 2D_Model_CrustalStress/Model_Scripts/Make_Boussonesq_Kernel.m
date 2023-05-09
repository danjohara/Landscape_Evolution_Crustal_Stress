function kernel = Make_Boussonesq_Kernel(y,z,stressComp,Fv,lambda,G,kernel_radius,kernel_res)
% Name: Make_Boussonesq_Kernel
% Python Code Authors: Richard H. Styron and Eric A. Hetland
% Matlab Adaptation: Daniel O'Hara
% Date: 11/20/2022
%
% Description: Matlab recreation of the Halfspace Python module from Styron
% & Hetland (2015)(https://github.com/cossatot/halfspace). Script generates
% the kernel used for the Boussonesq convolution. Used in a coupled 
% landscape evolution - crustal stress model by O'Hara & Karlstrom (in 
% review).
%
% Input:
%   y:              y-coordinate (currently always NaN for profiles).
%   z:              z-coordinate
%   stressComp:     Stress component to create kernel.
%   Fv:             Gravitational force.
%   lambda:         First Lame parameter.
%   G:              Second Lame parameter.
%   kernel_radius:  Kernel radius for convolution.
%   kernel_res:     Kernel resolution (grid resolution).
%
% Ouput:
%   kernel:         Stress kernel for convolution.

%% Make Kernel Array/Grid
kernel_len = round(kernel_radius*2/(kernel_res+1));

if isnan(y)
    kernelX = linspace(-kernel_radius,kernel_radius,kernel_len);
    kernelY = 0;
else
    tmpX = linspace(-kernel_radius,kernel_radius,kernel_len);
    [kernelX,kernelY] = meshgrid(tmpX,tmpX);
end

% Calculate Kernel
groundR = sqrt(kernelX.^2 + kernelY.^2);
R = sqrt(kernelX.^2 + kernelY.^2 + z.^2);
switch stressComp
    case 'XX'
        t1 = Fv./(2*pi);
        t2 = 3*kernelX.^2.*z./R.^5;
        t3 = G*(kernelY.^2 + z^2)./((lambda+G)*(z+R).*R.^3);
        t4 = (G*z)./((lambda+G).*R.^3);
        t5 = (G*kernelX.^2)./((lambda+G).*R.^2.*(z+R).^2);
        kernel = t1.*(t2 + t3 - t4 - t5);
    case 'YY'
        t1 = Fv./(2*pi);
        t2 = 3*kernelY.^2.*z./R.^5;
        t3 = G*(kernelX.^2 + z^2)./((lambda+G)*(z+R).*R.^3);
        t4 = (G*z)./((lambda+G).*R.^3);
        t5 = (G*kernelY.^2)./((lambda+G).*R.^2.*(z+R).^2);
        kernel = t1.*(t2 + t3 - t4 - t5);
    case 'ZZ'
        kernel = 3*Fv*z^3./(2*pi*R.^5);
    case 'XY'
        t1 = Fv./(2*pi);
        t2 = (3*kernelX.*kernelY.*z)./R.^5;
        t3 = G*kernelX.*kernelY.*(z+2*R)./((G+lambda)*R.^3.*(z+R).^2);
        kernel = t1.*(t2 - t3);
    case 'XZ'
        kernel = 3*Fv.*kernelX.*z^2./(2*pi*R.^5);
    case 'YZ'
        kernel = 3*Fv.*kernelY.*z^2./(2*pi*R.^5);
end

% Scale kernel for resolution
kernel = kernel*kernel_res^2;

% Make kernel circular if a grid.
if ~isnan(y)
    kernel(groundR > kernel_radius) = 0;
end
end